# Step 5 redo — eliminate per-protein DataFrames (numpy slices on the hot path)

This document distills PR #88 into its **final, intended design** and lays out a
clean, **step-by-step reviewable** implementation: a small number of commits that
are each logically complete, individually green (suite passes), and individually
bit-exact against the reference.

It supersedes the step-5 section of `IMPLEMENTATION_PLAN.md`, which describes an
earlier design that PR #88 reversed during review (it duplicated the normalization
math into `protein_intensity_estimation.py` and kept `ProtvalCutter`). This plan
bakes the final decisions in from the first commit, so the history reads as a
clean build instead of a build-then-unbuild.

> Steps 1–4, 6, 7 of the optimization series are unchanged — see
> `IMPLEMENTATION_PLAN.md`. Step 5 sits on top of step 4
> (`opt/step4-numpy-normalize-quadratic-linear`).

---

## 1. Why this step exists

After step 4 the per-protein numerical math is already numpy, but the **work items
crossing the multiprocessing boundary are still per-protein pandas DataFrames**:

- `get_normed_dfs` / `get_subdf` build ~14,240 `MultiIndex` sub-DataFrames
  **serially** (~2.5 s just to construct them).
- Each is **pickled to a worker and back** (~69 MB total); the equivalent numpy is
  ~10× cheaper to serialize. This is why multi-core scaling stalled at ~4× on 8 cores.

**Goal:** replace the per-protein pandas interface end to end with numpy. The
producer ships contiguous array slices, the worker computes on arrays and returns
arrays, and the assembly functions read those arrays.

**Measured gain:** estimate stage 16.2 s → **9.6 s** single-core (~1.7×), and 8-core
estimate ~10.5 s → **~2.4 s** (scaling restored; serial sub-DataFrame build 2.5 s → ~0.01 s).

---

## 2. Key design decisions (the distilled "what")

These are the decisions the final PR converged on. They are the contract this redo
must hit; the commit plan in §5 is just how to land them reviewably.

1. **Contiguous numpy slices, not DataFrames, are the unit of work.**
   `normed_df` is sorted by `(protein, ion)`, so each protein is a contiguous row
   block. `find_nameswitch_indices(protein_names)` gives the block boundaries; the
   producer slices `normed_df.to_numpy()` per protein. No per-protein DataFrame is
   ever built.

2. **The work item is a 6-tuple, the worker result a 4-tuple.**
   - Work item: `(idx, protein_name, ion_names, peptide_values, num_samples_quadratic, min_nonan)`
     — positional match to `calculate_peptide_and_protein_intensities`.
   - Result: `(protein_profile, protein_name, ion_names, shifted_values)`.
     `peptide_values`/`shifted_values` have **rows = ions, columns = samples**.
     `protein_profile` is `None` when every sample collapses to NaN.

3. **One shared normalization implementation, living in `normalization.py`.**
   The per-protein normalization is a module-level function
   `normalize_protein_ion_values(peptide_values, num_samples_quadratic)`.
   `NormalizationManagerProtein._normalize_quadratic_and_linear` *delegates* to it,
   and the hot path calls it directly. There is **no duplicated copy** in
   `protein_intensity_estimation.py`.
   - *Why this matters:* step 4 added a numpy override on the class; an early
     version of step 5 copied that logic into the hot path, leaving two divergent
     copies of order-sensitive float math. A single source removes that risk and is
     the reason `NormalizationManagerProtein` (still used by `visualizations.py`)
     stays correct.

4. **`normalize_protein_ion_values` reproduces *both* dispatch branches.**
   The base `NormalizationManager._run_normalization` dispatches on ion count:
   `n <= k` → quadratic-only (`_normalize_complete_input_quadratic`, rows in
   **original order**); `n > k` → quadratic+linear. The hot path bypasses that
   dispatch, so the standalone function must contain **both** branches. The
   class only ever reaches it on the `n > k` branch, so the class stays byte-identical.
   - The `n <= k` branch must **not** be folded into the sorted split — the
     quadratic-only path feeds rows to `get_normfacts` in original order, and
     `get_normfacts`' pair selection is order-sensitive.

5. **`ProtvalCutter` (class + numba `_get_num_nas_in_row`) is replaced and removed.**
   A numpy helper `_cut_peptide_values(peptide_values, ion_names, maximum=100)`
   does the same "keep ≤100 ions, sort by NaN-count asc then summed-intensity desc"
   reduction. `get_subdf` / `get_normed_dfs` are also removed. `OrphanIonRemover` /
   `OrphanIonsForDeletionSelector` / `IonCheckedForOrphan` (dead, only listed in
   `__all__`) go too.

6. **Single-sample guard is kept.** When `peptide_values.shape[1] < 2` (one sample),
   skip normalization and keep values as-is. Otherwise `get_normfacts` NaNs out
   every single-intensity ion row and the protein drops out of the output entirely.

---

## 3. The two bit-exactness gotchas (the correctness crux)

Float summation/median are **not associative**, so changing reduction *order* or
*memory layout* flips the last bit and propagates to a ~1e-6 table difference. Both
traps below were found by bisecting an `exact=False` verdict and must be reproduced
exactly by the numpy code.

**Gotcha 1 — `np.nansum` memory order.** The old `summed_pepint =
np.nansum(2 ** peptide_intensity_df)` summed a homogeneous float DataFrame, which
pandas stores **column-major (Fortran order)**. `np.nansum(2**c_contiguous_slice)`
sums row-major → differs by ~1 ULP, which flows through
`log2(summed_pepint / …)` into the **protein** table (the ion table never uses
`summed_pepint`). **Fix:** `np.nansum(np.asfortranarray(2 ** peptide_values))`.
Affects only proteins with >100 ions (12 on HYE), so unit data is too small to
catch it — the **optbench reference check is the real guard**.

**Gotcha 2 — `np.lexsort` must reproduce `sorted()`'s tie-break.** `ProtvalCutter`
sorted with `sorted(idxs, key=lambda i: (nan_count(i), -nansum(i)))` (stable).
`_cut_peptide_values` must produce the **same order**, not just the same set, or the
quadratic-subset selection *and* the `summed_pepint` order change. Use
`np.lexsort((neg_summed, nan_counts))` — last key is primary, lexsort is stable, so
full ties keep original (ion-name) order, matching `sorted()`.

---

## 4. Verification protocol (run after every commit)

1. **`python optbench.py`** → must print `VERDICT identical=True` with
   `max_abs_diff == 0` for both the protein and ion tables (sorted by id columns,
   compared with `equal_nan=True`). The baseline reference is captured once on the
   pre-step-5 code via `python optbench.py save`.
2. **pytest** (excluding seaborn-dependent viz/benchmark modules) — stays green.
3. Run from a neutral working directory, not the repo root (the `directlfq/` package
   dir would shadow the installed package).

Every commit below is expected to pass both 1 and 2.

---

## 5. Implementation plan — 5 reviewable commits

Each commit is logically complete, leaves the suite green, and stays bit-exact.
The ordering follows *"make the change easy, then make the easy change"*: a
behavior-preserving refactor (commit 1) first reshapes the code so the later numpy
flip is a small, obvious diff, then the helpers are proven against the code they
replace, then the interface flips, then the dead code is deleted.

> **Why the upfront refactor (commit 1) matters.** PR #88 changed two independent
> things in one commit: (a) how protein/ion **metadata reaches the assembly
> functions** — they re-derived `protein_name`/`ion_names` from the result
> DataFrame's `MultiIndex` — and (b) the **compute substrate** (DataFrame → numpy).
> Bundled, the diff is large and the assembly functions appear to change "because of
> numpy" when they really change because of metadata plumbing. Commit 1 does (a)
> alone, *while still on DataFrames*, so the assembly functions reach their final
> shape early and the numpy flip (commit 4) doesn't touch them at all.

### Commit 1 — Preparatory refactor: explicit metadata tuple (protein_intensity_estimation.py)

**What.** Behavior-preserving, still 100% DataFrame-based. Two changes:

1. **Rename** `get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan`
   → `get_protein_workitems` (pure rename; body unchanged, still builds DataFrames
   via `get_normed_dfs`/`get_subdf`).
2. **Make the worker return metadata explicitly.** Today
   `calculate_peptide_and_protein_intensities` returns `(protein_profile,
   shifted_df)`, and the assembly functions dig `protein_name`/`ion_names` out of
   `shifted_df.index`. Change the worker to extract them **after** the `ProtvalCutter`
   cut and return the final 4-tuple shape `(protein_profile, protein_name, ion_names,
   shifted_values)`, where `shifted_values = shifted_df.to_numpy()`. Update both
   assembly functions to read the tuple fields directly (`protein_name`,
   `ion_names.tolist()`, `np.concatenate(ion_vals)` over the arrays).

Everything numeric stays on DataFrames: the worker still uses `ProtvalCutter`,
`NormalizationManagerProtein`, and `np.nansum(2 ** peptide_intensity_df)` (no Fortran
fix yet — that arrives with the array switch in commit 4).

**Why it's safe.** No reduction order or memory layout changes; `shifted_df.to_numpy()`
carries the identical numbers the assembly already converted via `.to_numpy()`.
Bit-exact, suite green.

**Payoff.** After this commit the assembly functions are in their **final** form, so
commit 4 (the numpy flip) shrinks to: producer body, worker compute substrate, and
the two profile helpers — and the reviewer never has to disentangle metadata plumbing
from the numpy migration.

**Tests.** Update tests that unpack the worker result / drive the assembly functions
to the 4-tuple shape (this is the only test churn for the new tuple — commit 4 won't
revisit it).

**Verify.** optbench identical + pytest green.

---

### Commit 2 — Extract the shared `normalize_protein_ion_values` (normalization.py)

**What.** Lift the numpy body of `NormalizationManagerProtein._normalize_quadratic_and_linear`
(from step 4) into a module-level function
`normalize_protein_ion_values(peptide_values, num_samples_quadratic)`. The function
contains **both** dispatch branches (decision §2.4); the class method becomes a thin
delegator:

```python
def _normalize_quadratic_and_linear(self) -> None:
    df = self.complete_dataframe
    arr = normalize_protein_ion_values(
        df.to_numpy(dtype=float), self._num_samples_quadratic
    )
    self.complete_dataframe = pd.DataFrame(arr, index=df.index, columns=df.columns)
```

Add `normalize_protein_ion_values` to `__all__`. The input array is not mutated
(`.copy()` inside).

**Why this is safe.** The class only calls this method on the `n > k` branch, and
that branch is the unchanged step-4 code, so `NormalizationManagerProtein` output is
byte-for-byte identical. The `n <= k` branch is new code that only the (future) hot
path exercises. This is decision §2.3 + §2.4 landed up front, so it never has to be
"consolidated" later.

**Tests.** Pin `normalize_protein_ion_values` against the **base-class `.loc`
pipeline** it reimplements — drive a plain `NormalizationManager` with
`normalization_function = normalize_ion_profiles` and call `_run_normalization()`,
*not* the delegating subclass (otherwise the test is circular). Two cases:
- quadratic+linear (`n=5 > k=3`, mixed NaNs);
- quadratic-only (`n=3 <= k=5`; verifies the original-order branch).

**Verify.** optbench identical (class unchanged) + pytest green.

---

### Commit 3 — Add `_cut_peptide_values`, proven against the live `ProtvalCutter` (protein_intensity_estimation.py)

**What.** Add the numpy cutting helper (decision §2.5, Gotcha 2). Do **not** wire it
into the worker yet — `ProtvalCutter` is still present and still used.

```python
def _cut_peptide_values(peptide_values, ion_names, maximum=100):
    if peptide_values.shape[0] <= maximum:
        return peptide_values, ion_names
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        neg_summed = -np.nansum(peptide_values, axis=1)
    nan_counts = np.isnan(peptide_values).sum(axis=1)
    order = np.lexsort((neg_summed, nan_counts))[:maximum]  # Gotcha 2
    return peptide_values[order], ion_names[order]
```

**Why a separate commit.** This is the strangler pattern: introduce and *prove* the
replacement beside the original, then switch (commit 4), then delete (commit 5).
While both exist, the helper can be tested for **direct equivalence against the live
`ProtvalCutter`** on real inputs — the strongest possible check, and one that
disappears once the class is gone.

**Tests.** Equivalence vs. `ProtvalCutter` on (a) a `>maximum` case including a full
tie on both keys (asserts ion order, not just the set) and (b) a within-limit no-op.

**Verify.** optbench identical (helper unwired) + pytest green.

> Optional: commits 3 and 4 can be merged if a separate "unwired helper" commit is
> unwanted. Keeping them split gives the reviewer a focused diff for the subtle
> lexsort tie-break with a direct A/B test against the class it replaces.

---

### Commit 4 — Flip the hot path to numpy slices (protein_intensity_estimation.py)

**What.** The irreducible interface change. Because the contract changes from
DataFrame to numpy, **tests and implementation land in the same commit** (the new
shape can't be verified against the pre-change code).

- **Producer** — replace `get_input_specification_…` / `get_normed_dfs` / `get_subdf`
  with one generator `get_protein_workitems(normed_df, num_samples_quadratic,
  min_nonan)` yielding the 6-tuple over `find_nameswitch_indices` slices (decision
  §2.1, §2.2).
- **Worker** — `calculate_peptide_and_protein_intensities` takes the 6-tuple, returns
  the 4-tuple; uses `_cut_peptide_values` (when `>1` ion),
  `lfqnorm.normalize_protein_ion_values`, the single-sample guard (§2.6), and the
  Fortran-order `summed_pepint` (Gotcha 1).
- **Profile helpers** — `get_protein_profile_from_shifted_peptides` and
  `get_list_with_protein_value_for_each_sample` take a numpy array (drop the
  `df.to_numpy()` line).
- **Assembly** — *unchanged.* They already read the 4-tuple and operate on arrays
  after commit 1. This is the payoff of the upfront refactor.
- **Pool wrappers** — `get_list_with_sequential_processing` /
  `get_list_with_multiprocessing` only get their parameter renamed (generic
  `map`/`starmap`, no logic change).

`ProtvalCutter` / `get_subdf` / `get_normed_dfs` are now dead but **stay** this
commit — removing them is commit 5, keeping this diff scoped to the interface flip.
The commit-3 `_cut_peptide_values`-vs-`ProtvalCutter` tests still pass.

**Tests.** Update the `get_list_with_protein_value_for_each_sample` tests to pass
arrays. Add: single-sample values are kept (not NaN-ed) and the protein is retained.
(The worker-returns-4-tuple and ion-df-from-tuples tests already exist from commit 1.)

**Verify.** optbench `identical=True` (this is where Gotcha 1 actually bites —
unit data is too small) + pytest green.

---

### Commit 5 — Remove dead code, de-circularize tests, update probes

**What.**
- Delete `ProtvalCutter` (and its `@njit _get_num_nas_in_row`, and the
  now-unused `from numba import njit`), `get_subdf`, `get_normed_dfs`,
  `OrphanIonRemover`, `OrphanIonsForDeletionSelector`, `IonCheckedForOrphan`.
  Prune all of these from `__all__`; add `get_protein_workitems`.
- **De-circularize the commit-3 tests:** the `ProtvalCutter` reference is gone, so
  replace the equivalence assertions with **hardcoded expectations** (same inputs,
  expected ion order/values written out).
- **Benchmark probes:** move the old `ProtvalCutter` sort/cut into
  `benchmark_/optbench.py` as `protvalcutter_reference_sorted_index` /
  `protvalcutter_reference_cut`, and update `layoutprobe.py`, `orderprobe.py`,
  `sumprobe.py` to call those (they currently reach into the deleted class).

**Why last.** The deletions are only safe once nothing references the symbols
(verified: no callers in `directlfq/` or `tests/` after commit 4; the only
`NormalizationManagerProtein` consumer, `visualizations.py`, is untouched by step 5).

**Verify.** optbench identical + pytest green.

---

## 6. Commit summary

| # | file(s) | change | green & bit-exact |
|---|---------|--------|:--:|
| 1 | `protein_intensity_estimation.py` (+ test) | **prep refactor (DataFrame-based):** rename → `get_protein_workitems`; worker returns explicit `(profile, name, ion_names, shifted_values)`; assembly reads tuple fields | ✓ |
| 2 | `normalization.py` (+ test_normalization) | extract module-level `normalize_protein_ion_values` (both branches); class delegates | ✓ |
| 3 | `protein_intensity_estimation.py` (+ test) | add `_cut_peptide_values`, proven vs. live `ProtvalCutter` (unwired) | ✓ |
| 4 | `protein_intensity_estimation.py` (+ test) | flip hot path: numpy producer/worker + array profile helpers (assembly untouched), Gotcha 1 + single-sample guard | ✓ |
| 5 | `protein_intensity_estimation.py`, `optbench.py`, probes (+ test) | delete `ProtvalCutter`/`get_subdf`/`get_normed_dfs`/Orphan*, prune `__all__`, de-circularize tests, point probes at optbench reference | ✓ |

After commit 5: `optbench.py` prints `VERDICT identical=True`, pytest is green, and a
1/2/4/8-core `bench_dlfq.py` run confirms the restored multi-core scaling.
