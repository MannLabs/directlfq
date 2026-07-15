# directLFQ — memory & multi-core implementation plan

Self-contained plan to reproduce the **memory-reduction** (section A) and
**multi-core scaling** (section B) changes. Intended as the *only* input another
engineer/Claude instance needs to re-do them from the stated baseline.

These build on the earlier speed work (see `IMPLEMENTATION_PLAN.md`, which took the
`estimate_protein_intensities` stage from 137.6 s → ~1.6 s single-core by replacing
per-protein pandas with numpy and numba-compiling the clustering). That plan and
this one together describe the full optimization set.

---

## 0. Result summary

Test data: HYE1-6, generic wide TSV, 301,528 ions × 24 runs, 14,240 proteins.
Machine: 32 cores. Memory measured as **peak PSS** (proportional set size, summed
over the process tree — the true unique RAM; summed RSS over-counts copy-on-write
worker pages ~4-6× and must not be used).

| axis | before this plan | after |
|---|---|---|
| peak RAM, 1 core (default, bit-exact) | 0.74 GB | **0.60 GB** |
| peak RAM, 1 core (opt-in float32) | – | **0.47 GB** (−36% vs 0.74) |
| estimate stage, 8 cores | 1.75 s (slower than 1 core!) | **0.83 s** (scales, 2× vs 1c) |
| default `num_cores` behaviour | all cores (counterproductive) | sequential (opt-in parallel) |

Every default-path change is **bit-identical** to the pre-change output
(verified). float32 is opt-in and ~6e-6 relative-accuracy lossy.

---

## 1. Baseline / prerequisites

Assume the repo is at the post-speed-optimization state (in this project: commit
`6b7b92f`). Concretely, this plan assumes these already exist in
`directlfq/protein_intensity_estimation.py`:

- `find_nameswitch_indices(protein_names)` → array of per-protein block-start row
  indices (`switch`), length `n_proteins+1`. Each protein occupies contiguous rows
  `switch[i]:switch[i+1]` in the protein-sorted matrix.
- `_normalize_protein_values(peptide_values, num_samples_quadratic)` → numpy
  normalization of one protein's ion×sample matrix.
- `get_protein_profile_from_shifted_peptides(shifted_values, summed_pepints, min_nonan)`.
- `_nanmedian_axis0/1` (numba per-axis medians).
- `calculate_peptide_and_protein_intensities(idx, protein_name, ion_names, peptide_values, num_samples_quadratic, min_nonan)`
  returning `(profile, protein_name, ion_names, shifted_values)`, which internally
  cut >100-ion proteins, computed `summed_pepint = np.nansum(np.asfortranarray(2**peptide_values))`,
  normalized, and built the profile.
- Input is read via pyarrow in `directlfq/utils.py::_read_wide_tsv` (all-columns
  path) then `to_pandas()`.

The **memory high-water mark of the pipeline is the import stage** (read + sort +
dedup + log-transform), not estimate or write. On 1 core the estimate stage peaks
*lower* than import. This is why section A targets import first.

---

## 2. Verification protocol (run after every change)

1. **Bit-exact reference (once, from baseline code):** run the estimate stage on
   the dataset and save protein + ion outputs (full float64). A ready harness is
   `benchmarks/optbench.py` (`python optbench.py save`).
2. **After each default-path change:** `python optbench.py` → must print
   `VERDICT identical=True` (protein and ion `max_abs_diff == 0`, `equal_nan`).
   Also run the pytest suite (26 tests; exclude seaborn-dependent
   `test_visualizations`/`test_benchmarking`).
3. **Memory:** `benchmarks/memprofile.py <cores>` runs the full `run_lfq` in a
   child and samples the `/proc` process tree (RSS + PSS). Compare peak PSS.
4. **float32 accuracy:** compare to the float64 reference with a relative
   tolerance (expect max rel ≈ 6e-6), NOT bit-exact.
5. **Parallel bit-exactness:** run `estimate_protein_intensities(num_cores=4/8)`
   and compare to the reference — must be identical (order-independent, sort by id).

> Floating-point note: reductions are order-sensitive. `np.nansum` on a
> C-contiguous array differs at ~1 ULP from the same sum on an F-contiguous array;
> the existing `np.asfortranarray(2**peptide_values)` reproduces the historical
> column-major DataFrame sum and must be preserved.

---

## SECTION A — reduce absolute memory

### A1. Remove transient copies in the import stage  (bit-exact)
**File:** `directlfq/utils.py`. **Gain:** peak PSS 1c 0.74→0.60, 8c 1.10→0.97 GB.

**Motivation.** The import high-water mark comes from (a) holding the Arrow table
and the pandas copy simultaneously during the read, and (b) `index_and_log_transform`
chaining `.replace()` (one full copy) then `np.log2()` (a second full copy).

**A1a — free Arrow memory during conversion.** In `_read_wide_tsv`, the pyarrow
branch:
```python
# before
df = table.to_pandas()
# after
df = table.to_pandas(split_blocks=True, self_destruct=True)
del table
```

**A1b — in-place log transform.** Replace `index_and_log_transform_input_df`:
```python
# before
def index_and_log_transform_input_df(data_df):
    data_df = data_df.set_index([config.PROTEIN_ID, config.QUANT_ID])
    return np.log2(data_df.replace(0, np.nan))
# after
def index_and_log_transform_input_df(data_df):
    data_df = data_df.set_index([config.PROTEIN_ID, config.QUANT_ID])
    values = data_df.to_numpy(dtype=np.float64, copy=True)   # dtype becomes config-driven in A3
    values[values == 0] = np.nan
    np.log2(values, out=values)
    return pd.DataFrame(values, index=data_df.index, columns=data_df.columns)
```
**Exactness:** replace-0-with-NaN then log2 == set 0→NaN then log2 (log2(NaN)=NaN).

### A2. In-place `nan_to_num` for the ion table  (bit-exact)
**File:** `directlfq/protein_intensity_estimation.py::get_ion_intensity_dataframe_from_list_of_shifted_peptides`.
**Gain:** removes one ~57 MB transient copy during ion-table compilation (no change
to the global peak on this dataset; helps larger inputs).
```python
# before
merged_ions = 2**np.concatenate(ion_vals)
merged_ions = np.nan_to_num(merged_ions)
# after
merged_ions = 2**np.concatenate(ion_vals)
np.nan_to_num(merged_ions, copy=False)
```
**Rejected (do NOT do):** chunked `to_csv` for the ion table. Verified byte-identical
output but pandas already streams to the file handle → **no** measurable memory
change (write-region PSS 0.441↔0.443 GB). Not worth the complexity.

### A3. Opt-in float32 intensity path  (NOT bit-exact; default off)
**Files:** `directlfq/config.py`, `directlfq/utils.py`, `directlfq/lfq_manager.py`.
**Gain (with A4):** peak PSS 1c 0.60→0.47 GB. **Accuracy:** max rel ≈ 6e-6.

**Motivation.** LFQ is ratio-based; float64 is overkill. float32 halves every big
matrix. Must be opt-in to preserve the bit-exact default.

**A3a — config flag.** In `config.py`:
```python
import numpy as _np
INTENSITY_DTYPE = _np.float64
def set_intensity_dtype(use_float32: bool = False):
    global INTENSITY_DTYPE
    INTENSITY_DTYPE = _np.float32 if use_float32 else _np.float64
```
**A3b — single cast point.** In A1b's `index_and_log_transform_input_df`, change
`dtype=np.float64` → `dtype=config.INTENSITY_DTYPE`. Everything downstream
(per-protein slices, normalization, ion table) inherits the dtype. The per-protein
numba temporaries stay float64 (tiny; harmless upcast on assignment back to the
float32 arrays).
**A3c — wire the flag.** In `lfq_manager.run_lfq`, add `use_float32=False` to the
signature and call `config.set_intensity_dtype(use_float32=use_float32)` alongside
the other `config.set_*` calls.

Caveat:
  The algorithm leaves log space and does magnitude-mixing sum reductions to conserve total intensity, in protein_intensity_estimation.py:

  167:  summed_pepint    = np.nansum(np.asfortranarray(2**peptide_values))
  192:  summed_intensity = np.nansum(2**intens_vec)
  195:  intens_conversion_factor = summed_pepints / summed_intensity
  196:  scaled_vec = intens_vec + np.log2(intens_conversion_factor)

  If peptide_values is float32, then 2**peptide_values is float32 and np.nansum accumulates in a float32 accumulator. Linear intensities span many orders of magnitude within one protein, so summing e.g. 1e11 + 1e3 in float32 silently drops the small term (ULP at 1e11 is ~8000). That biases
  summed_pepint.

### A4. Read intensities directly as float32  (part of the opt-in path)
**File:** `directlfq/utils.py::_read_wide_tsv`. **Gain:** 1c peak 0.62→0.47 GB
(stable). Attacks the import peak, which A3 alone did not (it read float64 then cast).

**Motivation.** With A3 the CSV was still parsed to float64 then cast, so the
float64 matrix still materialized at the import peak. Parse the intensity columns
straight to float32.
```python
if samples_subset is None:
    try:
        import pyarrow as pa
        import pyarrow.csv as pacsv
        convert_options = None
        if config.INTENSITY_DTYPE == np.float32:
            with open(file_to_read) as fh:
                header = fh.readline().rstrip("\n").split("\t")
            col_types = {c: pa.float32() for c in header
                         if c not in (config.PROTEIN_ID, config.QUANT_ID)}
            col_types[config.PROTEIN_ID] = pa.string()
            col_types[config.QUANT_ID] = pa.string()
            convert_options = pacsv.ConvertOptions(column_types=col_types)
        table = pacsv.read_csv(file_to_read,
                               parse_options=pacsv.ParseOptions(delimiter="\t"),
                               convert_options=convert_options)
        df = table.to_pandas(split_blocks=True, self_destruct=True)
        del table
        return df
    except Exception as e:
        LOGGER.info(f"pyarrow CSV read failed ({e}); falling back to pandas read_csv.")
return pd.read_csv(file_to_read, sep="\t", encoding="latin1", usecols=samples_subset)
```
**Default path unchanged:** when `INTENSITY_DTYPE==float64`, `convert_options=None`
→ identical read, bit-exact.
**Note:** float32 helps the single-core (import) peak but **not** multi-core
(the fork/estimate transient there is dominated by numba code + COW of non-data
pages, and has ±0.1 GB run-to-run variance that swamps the float32 saving).

---
