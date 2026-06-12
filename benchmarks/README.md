# directLFQ speedup benchmarks

Scripts and notes from the single-core / multi-core optimization work on
`estimate_protein_intensities` and the sample-normalization step. See
**`NOTES.md`** for the full investigation, results, and commit-by-commit changes.

## Prerequisites

- The venv with `directlfq` editable-installed: `/home/magnus/dlfq_venv`.
- The converted HYE dataset at the absolute path hard-coded in the scripts:
  `/home/magnus/speedup_directlfq/directlfq/input_data/HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv`
  (produced by `hye_repro/convert_peaks_to_generic.py`).
- Run python from a neutral cwd (not the original repo dir, which would shadow
  the install). The probe scripts `import optbench`, so add this directory to
  `PYTHONPATH`:

```bash
cd /home/magnus
export PYTHONPATH=/home/magnus/dlfq_work/benchmarks
PY=/home/magnus/dlfq_venv/bin/python
```

## Main harnesses

| script | purpose |
|--------|---------|
| `bench_dlfq.py <cores>`   | stage-resolved full-pipeline timing (import / norm / estimate / wall) |
| `bench_estimate.py <cores>` | breaks the estimate stage into build / map / assemble / compile |
| `optbench.py [save]`      | times the single-core estimate **and** checks bit-exact identity of protein + ion outputs against a saved reference |

`optbench.py save` regenerates `opt_protein.reference.npz` / `opt_ion.reference.npz`
(kept next to the script). A bare `optbench.py` run compares the current code
against that reference. The reference shipped here was captured from the
**baseline** code, so a passing run proves the optimizations are bit-identical to
the original directLFQ.

```bash
$PY bench_dlfq.py 8
$PY optbench.py          # -> VERDICT identical=True
```

## One-off profiling probes

Each was used to pin down a specific bottleneck (see NOTES.md §4–§5b):

| script | what it measured |
|--------|------------------|
| `prof.py`        | cProfile of the sequential per-protein map (2500-protein subset) |
| `jit_probe.py`   | one-time numba JIT cost per process |
| `pickle_probe.py`| DataFrame vs numpy pickle cost (IPC) |
| `normprof.py`    | cProfile of the sample-normalization step |
| `diffprobe.py`   | which proteins differed after the numpy-slice refactor |
| `orderprobe.py`  | ProtvalCutter vs lexsort ordering |
| `sumprobe.py`    | `summed_pepint` old (DataFrame) vs new (array) |
| `layoutprobe.py` | found the C- vs Fortran-order `nansum` 1-ULP difference |

> Note: `opt_ion.reference.npz` is ~70 MB. If committing this directory to git,
> consider `.gitignore`-ing the `*.npz` reference files and regenerating them with
> `optbench.py save` instead.
