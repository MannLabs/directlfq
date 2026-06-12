# directLFQ performance optimization — implementation plan

This document is **self-contained**: it is intended to be the *only* input another
engineer (or Claude Code instance) needs to reproduce, from the unmodified
directLFQ source, the full set of performance optimizations described here.

It covers, for each step: the **motivation/reasoning**, the **exact functions and
files to change**, **before → after** code, the **correctness constraints**
(every change must be numerically bit-identical), and the **measured gain**.

---

## 0. TL;DR of the result

Optimizing the protein-intensity estimation and sample-normalization of directLFQ,
on the HYE1-6 benchmark (301,528 ions × 24 runs, 14,240 proteins):

| metric | before | after (1 core) | after (8 cores) |
|--------|-------:|---------------:|----------------:|
| `estimate_protein_intensities` stage | 137.6 s | **9.6 s** (14.3×) | 2.4 s |
| sample-normalization stage | 1.45 s | **0.75 s** (cached) | 0.75 s |
| full `run_directlfq.py` pipeline | 143.9 s | **14.2 s** (10.1×) | 7.4 s (19.4×) |

**Every change produces bit-identical protein + ion output** (`max_abs_diff == 0`,
`np.array_equal(..., equal_nan=True)`), and the repo's pytest suite (26 tests)
still passes. The recurring theme: the bottleneck was never the numerical math
(numba) — it was **per-protein pandas overhead** and **DataFrame serialization
across processes**. The fixes replace pandas with numpy on the hot path.

---

## 1. Background: what the code does and where the time goes

### Pipeline (`directlfq/lfq_manager.py::run_lfq`)
1. **import** — read the input table, sort by (protein, ion), drop duplicate ions,
   log2-transform, drop all-NaN rows. *Serial.*
2. **sample normalization** — `NormalizationManagerSamplesOnSelectedProteins`
   shifts the 24 samples onto a common scale. *Serial.*
3. **estimate** — `estimate_protein_intensities`: for each protein independently,
   normalize its ions and derive one protein-intensity profile. Parallelized
   per-protein over `multiprocess.Pool`. *~98% of original runtime.*
4. **write** result tables. *Serial.*

### Data layout (important for every change below)
After import, `normed_df` is a `pd.DataFrame` whose index is a **MultiIndex
(protein, ion)**, sorted by protein then ion, and whose 24 columns are the runs.
Because it is protein-sorted, each protein occupies a **contiguous block of rows**.
`find_nameswitch_indices(protein_names)` returns the block boundaries — this is the
key that lets us operate on plain numpy slices instead of per-protein DataFrames.

### Profiling findings (the motivation)
Profile the single-core estimate stage (`cProfile`, sort by `tottime`). On the
unmodified code the top entries are **all pandas internals**: `isinstance`,
`sanitize_array`, `Index.__new__`, `pandas...indexing._align_series`,
`MultiIndex.__getitem__`, `Series.__init__`. The numba kernels (`calc_nanmedian`,
`calc_nanvar`) do **not** appear near the top. Conclusions that drive the plan:

- Single-core cost is ~9.5 ms/protein on tiny 24-wide arrays → pandas wrapping
  overhead, not real math.
- Multiprocessing scaled only ~4× at 8 cores because **14,240 per-protein
  DataFrames are pickled to workers and back** (~69 MB; the equivalent numpy is
  ~10× cheaper to pickle) and the per-protein sub-DataFrames are **built serially**.
- A one-time **~0.6 s numba JIT** is paid in the (first) sample-normalization call.

---

## 2. Prerequisites & environment

1. **Install** directlfq in editable mode with the `stable` extra into a venv.
   `pip install -e ".[stable]"`. Pinned deps: numpy 1.23.5, pandas 2.2.3,
   numba 0.56.4, multiprocess 0.70.14, pyarrow 17.0.0.
   - `pyyaml` is imported by `directlfq/utils.py` but missing from the requirements
     files → `pip install pyyaml` (or add it to `requirements/requirements*.txt`).
   - `pytest`, `matplotlib`, `seaborn` are only needed to run the test suite.
2. **Dataset.** A generic wide-format directLFQ input (`*.aq_reformat.tsv`):
   columns `protein`, `ion`, then one intensity column per run. The reference
   numbers here use HYE1-6 (301,528 ions × 24 runs, 14,240 proteins).
3. **Run python from a neutral working directory**, not the repo root — the repo
   contains a `directlfq/` package dir that would shadow the installed package.

---

## 3. Verification protocol (do this FIRST, before any change)

All changes must be **bit-identical**. The protocol:

1. **Capture a reference from the UNMODIFIED code.** Run the estimate stage on the
   dataset once with the original source and save the protein + ion outputs.
2. After **each** change, re-run and compare to that reference: assert
   `np.array_equal(new, ref, equal_nan=True)` for both the protein table and the
   ion table (after sorting both by their id columns so row order is irrelevant).
3. Also run the pytest suite (excluding the seaborn-dependent viz/benchmarking
   modules) — it must stay green.

A ready-made harness is in this directory: **`optbench.py`**.
- `python optbench.py save` — generate the reference (run this on the *baseline*).
- `python optbench.py` — run + print `VERDICT identical=True/False` and
  `max_abs_diff` for protein and ion tables.

The essential comparison logic (in case you need to recreate it):
```python
# sort by id columns, then compare numeric matrices with equal_nan
def df_to_arrays(df):
    df = df.reset_index()
    id_cols  = [c for c in df.columns if c in ("protein", "ion")]
    val_cols = [c for c in df.columns if c not in id_cols and df[c].dtype.kind in "fiu"]
    ids  = df[id_cols].astype(str).agg("\t".join, axis=1).to_numpy()
    order = np.argsort(ids, kind="stable")
    return ids[order], df[val_cols].to_numpy(float)[order]
# exact = np.array_equal(vals, ref_vals, equal_nan=True)
```

> Why bit-exactness is non-trivial here: floating-point summation is **not
> associative**, so changing the *order* or *memory layout* of a reduction
> (`np.nansum`, `np.nanmedian`) can change the last bit. Two such traps are called
> out explicitly in steps 5 and the gotchas section.

---

## 4. The optimization steps

Apply in order. Files are under `directlfq/`. After each step, run the
verification protocol. Gains are the single-core `estimate` stage time on HYE1-6
(except #6/#7 which are the sample-normalization stage).

Add `import warnings` to the top of **both** `protein_intensity_estimation.py` and
`normalization.py` (used by several steps).

---

### Step 1 — Vectorize `get_list_with_protein_value_for_each_sample`
**File:** `directlfq/protein_intensity_estimation.py`
**Gain:** 137.6 → 123.9 s

**Motivation.** For every protein this builds one pandas `Series` per sample (24
per protein × 14,240 proteins) just to take a column median.

**Before:**
```python
def get_list_with_protein_value_for_each_sample(normalized_peptide_profile_df, min_nonan):
    intens_vec = []
    for sample in normalized_peptide_profile_df.columns:
        reps = normalized_peptide_profile_df.loc[:,sample].to_numpy()
        nonan_elems = sum(~np.isnan(reps))
        if(nonan_elems>=min_nonan):
            intens_vec.append(np.nanmedian(reps))
        else:
            intens_vec.append(np.nan)
    return intens_vec
```

**After (this step):** operate on `df.to_numpy()` in one pass — column-wise
nanmedian, then NaN out columns with `< min_nonan` finite ions.
```python
def get_list_with_protein_value_for_each_sample(normalized_peptide_profile_df, min_nonan):
    arr = normalized_peptide_profile_df.to_numpy()
    nonan_counts = np.sum(~np.isnan(arr), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN columns -> NaN
        intens_vec = np.nanmedian(arr, axis=0)
    intens_vec[nonan_counts < min_nonan] = np.nan
    return intens_vec
```
**Exactness.** `np.nanmedian(arr, axis=0)` computes exactly the per-column median
the loop did; masking reproduces the `nonan_elems >= min_nonan` branch. (This
function's signature changes again in step 5 to take a numpy array directly.)

---

### Step 2 — Vectorize `NormalizationManager._determine_sorted_rows`
**File:** `directlfq/normalization.py`
**Gain:** 123.9 → 116.9 s

**Motivation.** Sorts the ion rows by NaN-count using a per-row `.loc[idx]` lookup
inside `sorted()`'s key — O(n) pandas MultiIndex lookups per protein with >10 ions.

**Before:**
```python
def _determine_sorted_rows(self):
    rows = self.complete_dataframe.index
    self._rows_sorted_by_number_valid_values = sorted(
        rows, key=lambda idx: self._get_num_nas_in_row(self.complete_dataframe.loc[idx,:].to_numpy()))
```

**After:**
```python
def _determine_sorted_rows(self):
    rows = self.complete_dataframe.index
    nan_counts = np.isnan(self.complete_dataframe.to_numpy()).sum(axis=1)
    order = np.argsort(nan_counts, kind="stable")
    self._rows_sorted_by_number_valid_values = [rows[i] for i in order]
```
**Exactness.** `np.argsort(kind="stable")` preserves original order for equal
keys, matching Python's stable `sorted()`. Output is the same list of index labels.

---

### Step 3 — Vectorize `SampleShifterLinear` row shifting
**File:** `directlfq/normalization.py`
**Gain:** 116.9 → 102.8 s

**Motivation.** Shifts each row by a scalar via `self.ion_dataframe.iloc[row, :] +=
distance` — pandas Series alignment on every row (the `_align_series` hot spot,
~20% of estimate time).

**Before:**
```python
def _shift_columns_to_reference_sample(self):
    num_rows = self._ion_dataframe_values.shape[0]
    for row_idx in range(num_rows):
        self._shift_to_reference_sample(row_idx)

def _shift_to_reference_sample(self, row_idx):
    distance_to_reference = self._calc_distance(
        samples_1=self._reference_intensities, samples_2=self._ion_dataframe_values[row_idx,:])
    self.ion_dataframe.iloc[row_idx, :] += distance_to_reference
```

**After:** compute all row shifts at once and apply them with a single block write.
```python
def _shift_columns_to_reference_sample(self):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN rows -> NaN
        distances = np.nanmedian(self._reference_intensities - self._ion_dataframe_values, axis=1)
    shifted = self.ion_dataframe.to_numpy() + distances[:, None]
    self.ion_dataframe.iloc[:, :] = shifted
```
**Exactness.** `_calc_distance` was `nanmedian(reference - row)` (and `nan` if the
row is all-NaN — `np.nanmedian` already returns `nan` there). Each row's shift
depends only on its own original values + the reference, so order-independent;
the single block assignment yields identical values. (`_shift_to_reference_sample`
becomes unused; leave it or delete it.)

---

### Step 4 — Numpy `NormalizationManagerProtein._normalize_quadratic_and_linear`
**File:** `directlfq/normalization.py`
**Gain:** 102.8 → 16.2 s  ← **the single biggest win (6.3×)**

**Motivation.** For proteins with > `num_samples_quadratic` (=10) ions the base
class normalizes a "quadratic" subset and linearly shifts the rest, doing it via
**label-based** `self.complete_dataframe.loc[list_of_tuples, :] = df` assignments
on a MultiIndex — once per such protein. These round-trips (pandas
`_align_series` / `_setitem_single_column` / tuple handling) dominated runtime.

**What the base algorithm does (must be reproduced exactly):**
- sort rows by NaN-count (stable); first `k=num_samples_quadratic` = quadratic
  subset, the rest (in original order) = linear subset;
- normalize the quadratic subset with `normalize_ion_profiles` (= `get_normfacts`
  then `apply_sampleshifts`; `get_normfacts` first NaNs out rows with <2 finite
  values);
- reference sample = column-wise **median (skipna)** of the normalized subset;
- shift each linear row by `nanmedian(reference - row)`.

**After:** add this override to the `NormalizationManagerProtein` class (which
otherwise inherits `_run_normalization` from the base, so this method is only hit
for the > k case; the ≤ k case still uses the unchanged
`_normalize_complete_input_quadratic`):
```python
class NormalizationManagerProtein(NormalizationManager):
    def __init__(self, complete_dataframe, num_samples_quadratic):
        super().__init__(complete_dataframe, num_samples_quadratic)
        self.normalization_function = normalize_ion_profiles
        self._rows_sorted_by_number_valid_values = None
        self._run_normalization()

    def _normalize_quadratic_and_linear(self):
        df = self.complete_dataframe
        arr = df.to_numpy(dtype=float, copy=True)
        k = self._num_samples_quadratic

        nan_counts = np.isnan(arr).sum(axis=1)
        order = np.argsort(nan_counts, kind="stable")
        q_pos = order[:k]
        linear_mask = np.ones(arr.shape[0], dtype=bool)
        linear_mask[q_pos] = False
        lin_pos = np.flatnonzero(linear_mask)            # original order

        q_vals = arr[q_pos].copy()
        sample2shift = get_normfacts(q_vals)              # mutates single-intensity rows -> NaN
        q_normed = apply_sampleshifts(q_vals, sample2shift)
        arr[q_pos] = q_normed

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN slices -> NaN
            reference = np.nanmedian(q_normed, axis=0)
            if lin_pos.size:
                shifts = np.nanmedian(reference - arr[lin_pos], axis=1)
                arr[lin_pos] = arr[lin_pos] + shifts[:, None]

        self.complete_dataframe = pd.DataFrame(arr, index=df.index, columns=df.columns)
```
**Exactness.** Reuses the same `get_normfacts`/`apply_sampleshifts`; pandas
`.median(axis=0)` is skipna so it equals `np.nanmedian`; stable argsort matches the
base `_determine_sorted_rows`; the quadratic subset is kept in **sorted order**
(`order[:k]`, not re-sorted) because `get_normfacts`'s pair selection is
order-sensitive. Verified bit-identical.

> Note: `NormalizationManagerProtein` is also used by `visualizations.py`, so keep
> the class. Step 5 will stop the *hot path* from using it, but the override stays
> correct and live for other callers.

---

### Step 5 — Eliminate per-protein DataFrames (ship numpy slices, return numpy)
**File:** `directlfq/protein_intensity_estimation.py`
**Gain:** 16.2 → 9.6 s single-core, **and restores multi-core scaling** (estimate
at 8 cores: ~10.5 s → ~2.4 s; serial sub-DataFrame build 2.5 s → 0.01 s).

**Motivation.** Even after step 4, the work items were per-protein pandas
DataFrames: built serially (`get_normed_dfs`/`get_subdf`) and pickled to/from
workers. Replace the entire per-protein interface with numpy: ship contiguous
array slices, compute on arrays, return arrays.

This step rewrites several functions. The relevant downstream consumers are the
two assembly functions, which must be updated to read from the new tuple shape.

**5a. Work-item producer.** Yield lightweight tuples instead of DataFrames:
```python
def get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan(normed_df, num_samples_quadratic, min_nonan):
    protein_names = normed_df.index.get_level_values(0).to_numpy()
    ion_names = normed_df.index.get_level_values(1).to_numpy()
    normed_array = normed_df.to_numpy()
    switch = find_nameswitch_indices(protein_names)
    n_proteins = len(switch) - 1
    return zip(
        range(n_proteins),
        (protein_names[switch[i]] for i in range(n_proteins)),
        (ion_names[switch[i]:switch[i + 1]] for i in range(n_proteins)),
        (normed_array[switch[i]:switch[i + 1]] for i in range(n_proteins)),
        itertools.repeat(num_samples_quadratic),
        itertools.repeat(min_nonan),
    )
```
(`get_normed_dfs`/`get_subdf` become unused; leave them — they are in `__all__`.
`get_list_with_sequential_processing` / `get_list_with_multiprocessing` are generic
`map`/`starmap` wrappers and need no change.)

**5b. Worker.** New signature; operate on numpy; return a 4-tuple. Replaces the use
of `ProtvalCutter` and `NormalizationManagerProtein` on the hot path with numpy
helpers:
```python
def calculate_peptide_and_protein_intensities(idx, protein_name, ion_names, peptide_values, num_samples_quadratic, min_nonan):
    if peptide_values.shape[0] > 1:
        peptide_values, ion_names = _cut_peptide_values(peptide_values, ion_names, maximum=100)

    if config.LOG_PROCESSED_PROTEINS and (idx % config.LOG_PROCESSED_PROTEINS_INTERVAL == 0):
        LOGGER.info(f"lfq-object {idx}")
    # asfortranarray reproduces the column-major summation order of the old
    # np.nansum(2**peptide_intensity_df) path (pandas stores a homogeneous float
    # frame column-major) -> keeps summed_pepint bit-identical. SEE GOTCHA 1.
    summed_pepint = np.nansum(np.asfortranarray(2**peptide_values))

    shifted_values = _normalize_protein_values(peptide_values, num_samples_quadratic)
    protein_profile = get_protein_profile_from_shifted_peptides(shifted_values, summed_pepint, min_nonan)
    return protein_profile, protein_name, ion_names, shifted_values


def _cut_peptide_values(peptide_values, ion_names, maximum=100):
    """Numpy equivalent of ProtvalCutter: keep at most `maximum` ions, sorted by
    NaN count asc then summed intensity desc. Only reorders when > maximum ions."""
    if peptide_values.shape[0] <= maximum:
        return peptide_values, ion_names
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        neg_summed = -np.nansum(peptide_values, axis=1)
    nan_counts = np.isnan(peptide_values).sum(axis=1)
    order = np.lexsort((neg_summed, nan_counts))[:maximum]   # SEE GOTCHA 2
    return peptide_values[order], ion_names[order]


def _normalize_protein_values(peptide_values, num_samples_quadratic):
    """Numpy equivalent of NormalizationManagerProtein. Rows are ions, cols samples."""
    if peptide_values.shape[0] <= num_samples_quadratic:
        values = peptide_values.copy()
        sample2shift = lfqnorm.get_normfacts(values)
        return lfqnorm.apply_sampleshifts(values, sample2shift)

    arr = peptide_values.copy()
    k = num_samples_quadratic
    nan_counts = np.isnan(arr).sum(axis=1)
    order = np.argsort(nan_counts, kind="stable")
    q_pos = order[:k]
    linear_mask = np.ones(arr.shape[0], dtype=bool)
    linear_mask[q_pos] = False
    lin_pos = np.flatnonzero(linear_mask)

    q_vals = arr[q_pos].copy()
    sample2shift = lfqnorm.get_normfacts(q_vals)
    q_normed = lfqnorm.apply_sampleshifts(q_vals, sample2shift)
    arr[q_pos] = q_normed

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        reference = np.nanmedian(q_normed, axis=0)
        if lin_pos.size:
            shifts = np.nanmedian(reference - arr[lin_pos], axis=1)
            arr[lin_pos] = arr[lin_pos] + shifts[:, None]
    return arr
```
(`_normalize_protein_values` deliberately duplicates step-4's logic so the hot path
never touches pandas; `NormalizationManagerProtein` is left intact for other callers.)

**5c. Profile helper now takes an array** (was changed in step 1 to take a df):
```python
def get_protein_profile_from_shifted_peptides(shifted_values, summed_pepints, min_nonan):
    intens_vec = get_list_with_protein_value_for_each_sample(shifted_values, min_nonan)
    intens_vec = np.array(intens_vec)
    summed_intensity = np.nansum(2**intens_vec)
    if summed_intensity == 0:
        return None
    intens_conversion_factor = summed_pepints/summed_intensity
    return intens_vec + np.log2(intens_conversion_factor)

def get_list_with_protein_value_for_each_sample(shifted_values, min_nonan):
    nonan_counts = np.sum(~np.isnan(shifted_values), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        intens_vec = np.nanmedian(shifted_values, axis=0)
    intens_vec[nonan_counts < min_nonan] = np.nan
    return intens_vec
```

**5d. Assembly functions** — read from the 4-tuple `(profile, protein_name,
ion_names, shifted_values)` instead of a DataFrame index:
```python
def get_ion_intensity_dataframe_from_list_of_shifted_peptides(list_of_tuple_..., column_names):
    ion_names, ion_vals, protein_names = [], [], []
    for _, protein_name, inames, shifted_values in list_of_tuple_...:
        ion_names.extend(inames.tolist())
        ion_vals.append(shifted_values)
        protein_names.extend([protein_name] * len(inames))
    merged_ions = np.nan_to_num(2**np.concatenate(ion_vals))
    ion_df = pd.DataFrame(merged_ions)
    ion_df.columns = column_names
    ion_df["ion"] = ion_names
    ion_df["protein"] = protein_names
    return ion_df.set_index(["protein", "ion"])

# in get_protein_dataframe_from_list_of_protein_profiles, change only the allprots line:
#   allprots = [x[1].index.get_level_values(0)[0] for x in ...]   # OLD
    allprots = [x[1] for x in list_of_tuple_w_protein_profiles_and_shifted_peptides]  # NEW
```
**Exactness.** Ion table is built from the same shifted values and ion names (set +
values match regardless of within-protein order). Protein profiles need the
`asfortranarray` summation fix (Gotcha 1) and the lexsort-must-equal-`sorted()`
property (Gotcha 2). After both, output is bit-identical and 26 pytest tests pass.

---

### Step 6 — Vectorize `set_samples_with_only_single_intensity_to_nan`
**File:** `directlfq/normalization.py`
**Gain:** sample-normalization stage ~1.45 → ~1.32 s (also marginally helps estimate)

**Motivation.** Inside `get_normfacts`, this NaNs out any sample (row) with <2
finite values, using a per-row Python `sum(~isnan)` over the (very long) ion axis.

**Before:**
```python
def set_samples_with_only_single_intensity_to_nan(samples):
    for idx in range(len(samples)):
        sample = samples[idx]
        if sum(~np.isnan(sample)) <2:
            sample[:] = np.nan
```
**After:**
```python
def set_samples_with_only_single_intensity_to_nan(samples):
    counts = np.count_nonzero(~np.isnan(samples), axis=1)
    samples[counts < 2] = np.nan
```
**Exactness.** Same predicate (<2 finite per row), same in-place mutation. Shared
by sample-norm and per-protein estimate; both stay bit-identical.

---

### Step 7 — Enable numba `cache=True` on the normalization kernels
**Files:** `directlfq/normalization.py`, `directlfq/tracefilter.py`
**Gain:** sample-normalization stage ~1.32 → ~0.75 s on any run after the first
(persists the ~0.6 s one-time JIT compile across processes/runs).

**Motivation.** `calc_nanmedian` / `calc_nanvar` (used by `get_normfacts`) and
`check_connected_traces` are the first numba functions hit in a run. Their JIT
compile (~0.6 s) was paid on every process start.

**Change:** add `cache=True` to the decorators.
```python
# normalization.py
@njit(cache=True)
def calc_nanvar(fcdist):
    return np.nanvar(fcdist)

@njit(cache=True)
def calc_nanmedian(fcdist):
    return np.nanmedian(fcdist)

# tracefilter.py
@njit(cache=True)
def check_connected_traces(matrix, trace_idx, visited):
    ...
```
**Exactness.** Caching changes nothing numerically. Requires the package directory
(or `NUMBA_CACHE_DIR`) to be writable. The first run compiles + writes the cache;
subsequent runs load it.

> Caveat learned the hard way: `cache=True` is **not** worth it for the small
> per-protein estimate kernels (their JIT is ~0.15 s and is already amortized
> because pool workers are forked after the parent compiled). It *is* worth it for
> these normalization kernels because they are the first/large-vector compiles.

---

## 5. Critical correctness gotchas (read before steps 5)

These are floating-point–associativity traps. Both were found by getting
`exact=False` with a ~1e-6 difference and bisecting to the cause.

**Gotcha 1 — `np.nansum` memory order (column- vs row-major).**
The old code computed `summed_pepint = np.nansum(2 ** peptide_intensity_df)` on a
**DataFrame**. A homogeneous float DataFrame is stored **column-major**, so
`np.nansum` summed in Fortran order. A plain `np.nansum(2**array)` on a C-contiguous
numpy slice sums in row-major order → differs by ~1 ULP, which propagates through
`log2(summed_pepint/…)` to a ~1e-6 difference in the protein table (the ion table
is unaffected because it never uses `summed_pepint`). **Fix:**
`np.nansum(np.asfortranarray(2**peptide_values))`. Affects only the 12 proteins with
>100 ions on HYE, but the fix is universal and cheap.

**Gotcha 2 — `np.lexsort` must reproduce Python `sorted()`'s tie-breaking.**
`ProtvalCutter` sorted ion rows with `sorted(idxs, key=lambda i:(nan_count(i),
-nansum(i)))` (stable). The numpy replacement must produce the **same order**, not
just the same set, otherwise the quadratic-subset selection and the `summed_pepint`
sum order change. Use `np.lexsort((neg_summed, nan_counts))` — last key is primary,
and lexsort is stable so full ties keep original order, matching `sorted()`. Verify
on a >100-ion protein that the kept-ion order is identical.

**General rule:** preserve row/element order and the exact reduction functions
(`np.nanmedian`, `np.nansum`, skipna semantics). When in doubt, reuse the existing
numba/numpy primitives (`get_normfacts`, `apply_sampleshifts`) rather than
re-deriving the math.

---

## 6. Summary table (apply in this order)

| # | file · function(s) | change | gain (stage) |
|---|--------------------|--------|--------------|
| 1 | `protein_intensity_estimation.py` · `get_list_with_protein_value_for_each_sample` | per-sample Series loop → numpy column nanmedian | est 137.6→123.9 s |
| 2 | `normalization.py` · `NormalizationManager._determine_sorted_rows` | per-row `.loc` sort → numpy nan-count + stable argsort | est 123.9→116.9 s |
| 3 | `normalization.py` · `SampleShifterLinear._shift_columns_to_reference_sample` | per-row `.iloc +=` → one numpy block shift | est 116.9→102.8 s |
| 4 | `normalization.py` · `NormalizationManagerProtein._normalize_quadratic_and_linear` (new override) | drop label-based `.loc[tuples,:]=df` round-trips → numpy | est 102.8→16.2 s |
| 5 | `protein_intensity_estimation.py` · `get_input_specification_…`, `calculate_peptide_and_protein_intensities`, new `_cut_peptide_values`/`_normalize_protein_values`, `get_protein_profile_from_shifted_peptides`, `get_list_with_protein_value_for_each_sample`, `get_ion_intensity_dataframe_…`, `get_protein_dataframe_…` | eliminate per-protein DataFrames: numpy slices in, numpy out | est 16.2→9.6 s + multi-core restored |
| 6 | `normalization.py` · `set_samples_with_only_single_intensity_to_nan` | Python `sum()` loop → `np.count_nonzero(axis=1)` | norm 1.45→1.32 s |
| 7 | `normalization.py` (`calc_nanmedian`,`calc_nanvar`), `tracefilter.py` (`check_connected_traces`) | `@njit` → `@njit(cache=True)` | norm 1.32→0.75 s (cached) |

After all steps: verify `optbench.py` prints `VERDICT identical=True`, run the
pytest suite (excluding seaborn-dependent modules), and re-run a 1/2/4/8-core
scaling test (`bench_dlfq.py`) to confirm the wall time and restored scaling.

## 7. What was intentionally NOT changed (and why)
- The O(n²) sample-clustering math in `get_normfacts`/`calc_distance` — it is the
  genuine remaining cost, but numba-compiling the whole clustering core (Python
  dict bookkeeping, exact `argmin` tie-breaking, masked-nanmedian edge cases) is
  high-risk for bit-exactness and was deferred.
- A thread-pool + `nogil` backend — unnecessary once DataFrame pickling was removed.
- `ProtvalCutter` / `NormalizationManagerProtein` classes — kept intact (tests +
  `visualizations.py` depend on them); the hot path simply stops using them.
