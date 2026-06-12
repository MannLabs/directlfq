# directLFQ speedup — investigation notes

Date: 2026-06-12
Goal: explore options to speed up directLFQ; install `[stable]`, profile the HYE
example dataset, run a 1/2/4/8-core scaling experiment.

> Location note: the project repo `/home/magnus/speedup_directlfq/directlfq` and
> its `hye_repro/` are **read-only** (owned by `alphadia`), so this file, the
> venv, and the benchmark scripts live under `/home/magnus/` instead.

---

## 1. Environment & install

This machine is hostile to a normal install:

- System Python 3.10 has **no `pip`, no `ensurepip`**, and **no `uv`/`virtualenv`/`conda`**.
- **No passwordless `sudo`** (can't `apt install python3.10-venv`).
- The repo is **read-only**, so editable installs that write `*.egg-info` into it fail.

Working recipe:

```bash
# 1. venv without pip (plain `venv` fails on ensurepip)
python3 -m venv --without-pip /home/magnus/dlfq_venv

# 2. bootstrap pip from a downloaded file
#    (piping installers to `sh` is blocked by the sandbox classifier)
curl -LsS https://bootstrap.pypa.io/get-pip.py -o /home/magnus/get-pip.py
/home/magnus/dlfq_venv/bin/python /home/magnus/get-pip.py

# 3. writable working tree (repo is read-only; input_data isn't tracked so the
#    clone is small and keeps git history + the `speedup` branch for committing)
git config --global --add safe.directory /home/magnus/speedup_directlfq/directlfq
git config --global --add safe.directory /home/magnus/speedup_directlfq/directlfq/.git
git clone /home/magnus/speedup_directlfq/directlfq /home/magnus/dlfq_work
git config --global --add safe.directory /home/magnus/dlfq_work
/home/magnus/dlfq_venv/bin/python -m pip install -e "/home/magnus/dlfq_work[stable]"

# 4. missing dependencies (pyyaml is a real gap, see bug below; the others are
#    only needed to run the pytest suite, not directlfq itself)
/home/magnus/dlfq_venv/bin/python -m pip install pyyaml
/home/magnus/dlfq_venv/bin/python -m pip install pytest matplotlib   # for tests/
```

Note: run python from a neutral cwd (e.g. `/home/magnus`), not from inside the
original repo dir — that dir contains a `directlfq/` package that would shadow the
editable install.

`[stable]` pins resolved exactly: numpy 1.23.5, pandas 2.2.3, numba 0.56.4,
multiprocess 0.70.14, pyarrow 17.0.0. directlfq 0.3.4-dev0.

### ❗ Bug: PyYAML missing from requirements

`directlfq/utils.py` does `import yaml`, but `pyyaml` is in **neither**
`requirements/requirements.txt` (stable) **nor** `requirements_loose.txt`. A clean
`pip install directlfq[stable]` therefore crashes on first import:

```
ModuleNotFoundError: No module named 'yaml'
```

Fix: add `pyyaml` to both requirements files. One-line PR, correctness not speed.

---

## 2. Dataset & pipeline

- Input: PEAKS DIA export `input_data/HYE1-6_repeat/lfq.dia.features.csv` (242 MB).
- `hye_repro/convert_peaks_to_generic.py` reshapes it into directLFQ's generic wide
  format (`protein`, `ion`, one intensity column per run) saved as
  `lfq.dia.features.csv.aq_reformat.tsv`.
- Dimensions: **301,528 ions × 24 runs, 14,240 proteins**.
- Run via `hye_repro/run_directlfq.py <cores>` → `lfq_manager.run_lfq(...)`.

### Pipeline stages (`lfq_manager.run_lfq`)

1. `import_data` + sort + dedup + log-transform + drop-all-NaN — **serial**
2. between-sample normalization (`NormalizationManagerSamplesOnSelectedProteins`) — **serial**
3. `estimate_protein_intensities` — **parallel** (per-protein over `multiprocess.Pool`)
4. write result tables — serial

---

## 3. Scaling experiment (1/2/4/8 cores)

Stage-resolved, **fresh process per core count** (each pays its own JIT — matches
real CLI usage). Harness: `/home/magnus/bench_dlfq.py`. Machine has 32 cores.

| cores | import (s) | norm (s) | estimate (s) | wall (s) | wall speedup | parallel eff. |
|------:|----------:|---------:|-------------:|---------:|-------------:|--------------:|
| 1     | 0.9       | 1.5      | 138.3        | 140.7    | 1.00×        | 100%          |
| 2     | 0.9       | 1.4      | 81.2         | 83.6     | 1.68×        | 84%           |
| 4     | 1.0       | 1.4      | 47.4         | 49.8     | 2.83×        | 71%           |
| 8     | 0.9       | 1.4      | 32.9         | 35.2     | 3.99×        | 50%           |

(Plain `run_directlfq.py 1` with default per-protein logging on: 143.9s — consistent.)

**`estimate_protein_intensities` is ~98% of runtime.** Scaling is clearly sublinear:
8 cores buys only ~4×.

---

## 4. Where the time goes / why it's sublinear

Import+norm is only ~2.4s, so **Amdahl's serial fraction is not the explanation**.
The loss is *inside* the estimate stage. Breakdown (`/home/magnus/bench_estimate.py`):

| sub-step                         | 1 core (s) | 8 cores (s) |
|----------------------------------|-----------:|------------:|
| build per-protein sub-DataFrames | 2.54       | 2.53        |  ← serial
| map (the actual compute)         | 135.4      | 24.6        |  ← parallel
| assemble protein dataframe       | 0.13       | 0.51        |  ← serial
| compile ion dataframe            | 0.51       | 0.56        |  ← serial

Probes (`/tmp/jit_probe.py`, `/tmp/pickle_probe.py`):

- **numba JIT recompile per worker ≈ 0.15 s/process → negligible** *for the
  estimate-stage per-protein kernels* (small functions). NOTE: this is narrower
  than it first looked — the *sample-normalization* kernels
  (`calc_nanmedian`/`calc_nanvar` over long vectors) cost ~0.6 s to compile, and
  `@njit(cache=True)` there *is* worth it (see change #7). Ruled out only for the
  per-worker estimate path.
- The map itself scales 135.4→24.6 s = **5.5× at 8 cores** (69% efficiency). The
  remaining loss is **IPC/pickling + load imbalance**, not serial Python.
- Each protein is wrapped in a pandas DataFrame with a MultiIndex
  (`get_subdf`/`get_normed_dfs`) and pickled to/from workers:
  **14,240 DataFrames = 69 MB, ~0.56 s to pickle + 1.11 s to unpickle**;
  the equivalent raw numpy is **58 MB and ~0.06 s** (≈10× cheaper).
- Single-core cost is **138 s / 14,240 ≈ 9.5 ms per protein** on tiny 24-wide
  arrays — that is pandas-wrapping overhead, not real math.

**Key takeaway:** the dominant cost is per-protein pandas overhead in a 14k-iteration
loop, and the parallel inefficiency is serialization of pandas objects. Reducing the
per-protein cost helps *both* the single-core time and the multi-core efficiency.

---

## 5. Options to speed it up (ranked by leverage / effort)

> What was actually implemented is in §5b. Option 1 (remove pandas) and option 3
> (cut serialization) were done across changes #1–#5; the `pyyaml` fix (option 4)
> is noted but left as a suggested PR; option 2 (thread pool + `nogil`) was **not**
> needed once the pandas overhead and DataFrame pickling were gone.

1. **Remove per-protein pandas overhead (highest leverage).** The hot loop allocates
   a MultiIndex DataFrame per protein and runs pandas ops on ~24-element arrays. The
   data is already a contiguous, protein-sorted numpy array with known name-switch
   indices (`find_nameswitch_indices`), so the kernel can operate on plain numpy
   slices. Cuts the single-core cost *and* the pickle payload. Likely beats adding cores.

2. **Thread pool + `@njit(nogil=True)` instead of process pool.** If the heavy kernels
   release the GIL, threads avoid all inter-process pickling and the serial
   sub-DataFrame build — directly attacks the 50% efficiency loss at 8 cores.
   Requires the per-protein work to be GIL-free (depends on #1 to strip pandas).

3. **Cut serialization without a full rewrite.** Ship `(start, end)` index ranges into
   one shared array rather than 14,240 pickled DataFrames; return numpy, assemble once.
   Also pass an explicit `chunksize` to `starmap`.

4. **Fix the `pyyaml` dependency gap** (correctness; trivial).

5. `@njit(cache=True)` — ~0.15 s for the per-protein estimate kernels (not worth
   it there), but ~0.6 s for the sample-normalization kernels — see change #7,
   where it *was* applied.

---

## 5b. Implemented optimizations & results

Working tree: **`/home/magnus/dlfq_work`** (clone of the read-only repo, branch
`speedup`; editable-installed into the venv). Harness: `/home/magnus/optbench.py`
(`save` writes a reference from baseline code; bare run times + checks bit-exact
identity of protein **and** ion outputs). Reference captured from baseline code,
so every later run verifies prep + estimate against the original.

All seven changes keep outputs **bit-identical** (`exact=True`, max_abs_diff 0 on
both protein and ion tables, 14,235 proteins / 294,778 ions; 26 pytest tests pass).
The recurring fix was replacing **pandas** overhead with numpy — the profile showed
the numba math was never the bottleneck — plus eliminating per-protein DataFrame
serialization (#5) and the one-time numba JIT (#7).

| # | change | estimate-stage time | commit |
|---|--------|--------------------:|--------|
| 0 | baseline | 137.6s | — |
| 1 | vectorize `get_list_with_protein_value_for_each_sample` (drop 24 Series/protein) | 123.9s | abea967 |
| 2 | vectorize `NormalizationManager._determine_sorted_rows` (numpy nan-count + stable argsort) | 116.9s | 157edd2 |
| 3 | vectorize `SampleShifterLinear` (kill per-row `.iloc[r,:] +=` align) | 102.8s | 628184a |
| 4 | numpy `NormalizationManagerProtein._normalize_quadratic_and_linear` (no label `.loc[tuples,:]=df` round-trips) | 16.2s | 515a7c1 |
| 5 | eliminate per-protein DataFrames: ship numpy slices to workers, return numpy, numpy `_cut_peptide_values`/`_normalize_protein_values` | **9.6s** | 26eea9c |
| 6 | vectorize `set_samples_with_only_single_intensity_to_nan` (sample-norm) | (norm 1.45→1.32s) | 61ded54 |
| 7 | numba `cache=True` on `calc_nanmedian`/`calc_nanvar`/`check_connected_traces` | (norm 1.32→0.75s cached) | 9ac3a16 |

- **Estimate stage: 137.6s → 9.6s (14.3×) single-core.** Change #4 alone was 6.3×
  — the label-based `.loc[list_of_tuples, :] = df` assignments on a MultiIndex,
  once per >10-ion protein, were the dominant cost.
- **Full CLI pipeline (`run_directlfq.py 1`): 143.9s → 14.2s (10.1×)** after all of
  #1–#7 (15.3s after #1–#5, before the sample-norm fixes). 8-core: 35.2s → 7.4s.
  The figures in change #6/#7 below isolate the sample-norm step itself.
- Change #5 bit-exactness note: `summed_pepint` had to use
  `np.nansum(np.asfortranarray(2**values))` — the old `np.nansum(2**DataFrame)`
  summed column-major (pandas stores a homogeneous float frame F-contiguous), so a
  plain C-order numpy sum differed by 1 ULP on the 12 >100-ion proteins. With the
  fortran-order sum, protein + ion outputs are bit-identical and 26 pytest tests pass.

### Re-scaling after the IPC/serial fixes (#5)

Eliminating the serial `build_subdfs` (2.5s → 0.01s) and shipping numpy instead of
pickled DataFrames restored multi-core scaling (which had collapsed to ~1.5× once
the per-protein pandas work was gone):

| cores | estimate (s) | wall (s)\* | estimate speedup | wall vs orig 1-core |
|------:|-------------:|-----------:|-----------------:|--------------------:|
| 1 | 9.57 | 11.99 | 1.00× | 11.7× |
| 2 | 6.48 | 8.93 | 1.48× | 15.8× |
| 4 | 3.66 | 6.00 | 2.61× | 23.5× |
| 8 | 2.48 | 4.93 | 3.86× | 28.5× |

\*wall excludes result-writing; full `run_directlfq.py`: 1-core 15.3s, 8-core 8.15s.
estimate-stage map now scales 4.5× at 8 cores (was 1.9×).

### Sample-normalization (serial prep) — changes #6, #7

The serial prep (~1.45s) was: ~0.6s one-time numba JIT (these kernels are the
first numba calls in a run), ~0.25s in a pure-Python `sum()` masking loop, and
~0.6s genuine O(24²) median/variance clustering over the long ion vectors.
#6 vectorized the masking loop; #7 added `cache=True` so the JIT is paid only on
the first-ever run. After a cache-warming run, norm is **0.75s** (from 1.45s).
The median/variance math is left as-is (algorithmic).

### Scaling after #6/#7 (cached numba)

| cores | import | norm | estimate | wall\* | full CLI (incl. write) |
|------:|-------:|-----:|---------:|-------:|-----------------------:|
| 1 | 0.89 | 0.75 | 9.28 | 10.92 | 14.18s |
| 2 | 0.89 | 0.75 | 6.31 | 7.95 | — |
| 4 | 0.95 | 0.74 | 3.60 | 5.29 | — |
| 8 | 0.89 | 0.74 | 2.37 | 4.01 | 7.40s |

### What's left (not done — higher risk / lower ROI)

1. The per-protein **O(n²) ion-clustering math** (`get_normfacts`/`calc_distance`)
   is now the bulk of both the estimate map and the residual ~0.6s of sample-norm.
   Squeezing it means numba-compiling the whole clustering core (Python-dict
   bookkeeping, exact argmin tie-breaking, masked-nanmedian edge cases preserved) —
   substantial and risky.
2. Residual IPC: results still pickle back through the pool; a fork-inherited
   shared-memory input array would remove the input-side copy entirely.
3. The remaining ~0.9s **import** (read TSV + sort/dedup/log-transform) is now the
   single largest serial item at high core counts.

---

## 6. Reproduction artifacts

All benchmark scripts now live in this directory (`dlfq_work/benchmarks/`); see
`README.md` here for how to run them (venv path, dataset path, `PYTHONPATH`).

- venv: `/home/magnus/dlfq_venv`
- writable git working tree (editable install target, branch `speedup` w/ opt commits): `/home/magnus/dlfq_work`
- `optbench.py [save]` — single-core estimate timing + bit-exact correctness check vs reference (`opt_*.reference.npz`, kept here)
- `bench_dlfq.py <cores>` — stage-resolved full-pipeline timing
- `bench_estimate.py <cores>` — estimate-stage sub-step breakdown
- `prof.py` — cProfile of the sequential per-protein map (2500-protein subset)
- `normprof.py` — cProfile of the sample-normalization step
- `jit_probe.py` — one-time numba JIT cost per process
- `pickle_probe.py` — DataFrame vs numpy pickle cost
- `diffprobe.py` / `orderprobe.py` / `sumprobe.py` / `layoutprobe.py` — pinned the bit-exactness issues during the numpy-slice refactor

Converted input (writable): `input_data/HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv`

Validation: `hye_repro/compare_outputs.py <new> <ref>` checks protein-id sets,
NaN patterns, and max abs/rel diff against a reference table — use it to confirm any
optimization keeps results identical.
