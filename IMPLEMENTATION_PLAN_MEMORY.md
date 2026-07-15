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

## SECTION B — reduce multi-core scaling overhead

### B.0 Why (measured)
After the speed work, per-protein compute is ~0.1 ms, so `multiprocess.Pool`
overhead dominates. Benchmarked (map-only wall, `benchmarks/explore_parallel.py`,
all bit-identical to sequential):

| strategy | 1c | 2c | 4c | 8c |
|---|--:|--:|--:|--:|
| sequential | 1263 ms | – | – | – |
| Pool.starmap (old) | 3196 | 1967 | 1641 | 1568 |
| chunksize=1 | 11963 | 5692 | 4642 | 7519 |
| **shared-memory + ranges** | 2885 | 1570 | 924 | **571** |
| ranges, no ion return | 1665 | 930 | 516 | 335 |
| threads (nogil kernels) | 1270 | 1045 | 1216 | 2953 |

Conclusions: (1) the old pool is **slower than sequential at every core count** —
it pickles ~58 MB of value slices out and ~58 MB of results back; (2) chunksize is
not the lever (default already auto-chunks; forcing 1 is catastrophic); (3)
**threads are a dead end** — the per-protein Python/numpy glue holds the GIL, so
even `nogil` kernels only help at 2c then degrade; (4) **shared-memory + index
ranges wins** (2.2× vs sequential at 8c) and is portable (works on fork and spawn).

### B1. Implement shared-memory + index ranges  (bit-exact)
**File:** `directlfq/protein_intensity_estimation.py`. **Gain:** estimate stage 8c
1.75→0.83 s; multi-core now *scales* instead of regressing.

**Design:** intensity matrix → a `multiprocessing.shared_memory` segment (workers
attach, not pickle); workers get only `(start,end)` protein-index ranges + the
small `switch` array via the pool initializer; protein/ion names stay in the
parent and are attached after workers return (no string metadata shipped).

**B1a — imports:** add `from multiprocessing import shared_memory as _shared_memory`.

**B1b — extract a names-independent core** so both sequential and parallel paths
share it, and expose the >100-ion cut permutation so the parent can reorder names:
```python
def calculate_peptide_and_protein_intensities(idx, protein_name, ion_names, peptide_values, num_samples_quadratic, min_nonan):
    if config.LOG_PROCESSED_PROTEINS and (idx % config.LOG_PROCESSED_PROTEINS_INTERVAL == 0):
        LOGGER.info(f"lfq-object {idx}")
    protein_profile, shifted_values, cut_order = _compute_protein_core(peptide_values, num_samples_quadratic, min_nonan)
    if cut_order is not None:
        ion_names = ion_names[cut_order]
    return protein_profile, protein_name, ion_names, shifted_values

def _compute_protein_core(peptide_values, num_samples_quadratic, min_nonan):
    cut_order = _cut_order(peptide_values, maximum=100)
    if cut_order is not None:
        peptide_values = peptide_values[cut_order]
    summed_pepint = np.nansum(np.asfortranarray(2**peptide_values))   # keep F-order (bit-exact)
    shifted_values = _normalize_protein_values(peptide_values, num_samples_quadratic)
    protein_profile = get_protein_profile_from_shifted_peptides(shifted_values, summed_pepint, min_nonan)
    return protein_profile, shifted_values, cut_order

def _cut_order(peptide_values, maximum=100):
    if peptide_values.shape[0] <= maximum:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        neg_summed = -np.nansum(peptide_values, axis=1)
    nan_counts = np.isnan(peptide_values).sum(axis=1)
    return np.lexsort((neg_summed, nan_counts))[:maximum]   # stable, matches ProtvalCutter order
```
(Keep a `_cut_peptide_values(values, names, maximum)` helper that applies
`_cut_order` to both, for backward compatibility.)

**B1c — worker + pool plumbing** (all module-level so they are picklable / present
after fork or spawn):
```python
def _split_ranges(n_items, n_chunks):
    edges = np.linspace(0, n_items, n_chunks + 1).astype(int)
    return [(int(edges[k]), int(edges[k+1])) for k in range(n_chunks) if edges[k] < edges[k+1]]

_PARALLEL_STATE = {}

def _shm_worker_init(shm_name, shape, dtype, switch, num_samples_quadratic, min_nonan):
    shm = _shared_memory.SharedMemory(name=shm_name)
    _PARALLEL_STATE.update(shm=shm,   # keep ref so the segment stays mapped
        arr=np.ndarray(shape, dtype=dtype, buffer=shm.buf),
        switch=switch, num_samples_quadratic=num_samples_quadratic, min_nonan=min_nonan)

def _shm_worker_range(bounds):
    lo, hi = bounds
    st = _PARALLEL_STATE
    arr, switch = st["arr"], st["switch"]
    nsq, min_nonan = st["num_samples_quadratic"], st["min_nonan"]
    return [_compute_protein_core(arr[switch[i]:switch[i+1]], nsq, min_nonan) for i in range(lo, hi)]

def get_list_with_shared_memory(normed_array, switch, protein_names, ion_names, num_samples_quadratic, min_nonan, num_cores):
    n_proteins = len(switch) - 1
    shm = _shared_memory.SharedMemory(create=True, size=max(int(normed_array.nbytes), 1))
    try:
        buf = np.ndarray(normed_array.shape, dtype=normed_array.dtype, buffer=shm.buf)
        buf[:] = normed_array
        bounds = _split_ranges(n_proteins, num_cores * 4)   # 4 chunks/core for load balance
        pool = multiprocess.Pool(num_cores, initializer=_shm_worker_init,
            initargs=(shm.name, normed_array.shape, normed_array.dtype, switch, num_samples_quadratic, min_nonan))
        try:
            parts = pool.map(_shm_worker_range, bounds)
        finally:
            pool.close(); pool.join()
    finally:
        shm.close(); shm.unlink()
    results, idx = [], 0
    for part in parts:
        for protein_profile, shifted_values, cut_order in part:
            start = switch[idx]                         # <-- CRITICAL: index by switch, see gotcha
            inames = ion_names[start:switch[idx+1]]
            if cut_order is not None:
                inames = inames[cut_order]
            results.append((protein_profile, protein_names[start], inames, shifted_values))
            idx += 1
    return results
```

**B1d — route the top-level dispatcher:**
```python
def get_list_of_tuple_w_protein_profiles_and_shifted_peptides(normed_df, num_samples_quadratic, min_nonan, num_cores):
    num_cores = _resolve_num_cores(num_cores)
    if num_cores <= 1:
        spec = get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan(normed_df, num_samples_quadratic, min_nonan)
        return get_list_with_sequential_processing(spec)
    protein_names = normed_df.index.get_level_values(0).to_numpy()
    ion_names = normed_df.index.get_level_values(1).to_numpy()
    normed_array = normed_df.to_numpy()
    switch = find_nameswitch_indices(protein_names)
    return get_list_with_shared_memory(normed_array, switch, protein_names, ion_names,
                                       num_samples_quadratic, min_nonan, num_cores)
```

> **CRITICAL GOTCHA (cost hours if missed).** `protein_names` and `ion_names` are
> per-**row** arrays (length = n_ions = 301,528), NOT per-protein. Protein `idx`'s
> name is `protein_names[switch[idx]]` (its first row), exactly as the sequential
> spec does with `protein_names[switch[i]]`. Using `protein_names[idx]` silently
> corrupts both the protein names and the ion table's protein index (it looked
> "protein max_abs_diff=0" but nan-pattern differed and the ion table was fully
> wrong). Always index names by `switch[idx]`.

### B2. Fix the default  (behaviour change)
**File:** same. Add and use `_resolve_num_cores`:
```python
def _resolve_num_cores(num_cores):
    if num_cores is None:
        return 1     # default = sequential; multi-core is opt-in (num_cores>=2)
    return num_cores
```
**Motivation.** The old default (`None` → all cores via `cpu_count`) was
counterproductive: the pool never beat sequential, and it used more memory. Per-
protein work is too cheap to parallelize by default. Sequential is fastest and
lightest for typical sizes; users with much larger inputs opt into `num_cores>=2`
and now get a real, bit-exact speedup via the shared-memory path.

**Not done (lower value):** a shared-memory *output* buffer to also eliminate the
result-side pickle (~40% of the residual IPC — see the `ranges, no ion return`
row). Results currently still pickle back to the parent.

---

## 3. Order, gains, and exactness — summary

| step | file · method(s) | gain | exact? |
|---|---|---|---|
| A1 | utils: `_read_wide_tsv`, `index_and_log_transform_input_df` | peak 1c 0.74→0.60 GB | yes |
| A2 | protein_intensity: `get_ion_intensity_dataframe_from_list_of_shifted_peptides` | −1 ion-matrix copy | yes |
| A3 | config: `INTENSITY_DTYPE`/`set_intensity_dtype`; utils cast point; lfq_manager `run_lfq(use_float32=)` | (enables A4) | default yes / float32 ~6e-6 |
| A4 | utils: `_read_wide_tsv` float32 column types | peak 1c 0.60→0.47 GB | default yes / float32 ~6e-6 |
| B1 | protein_intensity: `_compute_protein_core`, `_cut_order`, `_split_ranges`, `_shm_worker_init`, `_shm_worker_range`, `get_list_with_shared_memory`, `get_list_of_tuple_…`, import shared_memory | estimate 8c 1.75→0.83 s | yes (seq + 2/4/8c) |
| B2 | protein_intensity: `_resolve_num_cores` | default no longer regresses | n/a (default→sequential) |

**Rejected options (documented so they are not re-attempted):** chunked `to_csv`
(no memory benefit; pandas already streams); threads / `nogil` kernels (GIL-bound
per-protein glue; do not scale); larger `starmap` chunksize (not the bottleneck —
data volume is); a fully-numba-batched estimate that would fold `summed_pepint`
into numba (`numba.nansum` ≠ numpy F-order `nansum`, breaks bit-exactness).

## 4. Reproduction harnesses (in `benchmarks/`)
- `optbench.py [save]` — estimate timing + bit-exact check vs reference.
- `memprofile.py <cores>` — full-pipeline peak RSS/PSS via `/proc` sampling
  (set env `DLFQ_F32=1` to exercise float32).
- `explore_parallel.py` — times all multi-core strategies, checks each bit-exact.
- Companion result write-ups: `MEMORY_PROFILE.md`, `PARALLEL_EXPLORATION.md`.
