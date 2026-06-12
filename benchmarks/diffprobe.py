import numpy as np, optbench
import directlfq.protein_intensity_estimation as lfqprot
normed = optbench.prep()
prot, ion = lfqprot.estimate_protein_intensities(normed, min_nonan=1, num_samples_quadratic=10, num_cores=1)
ids, vals = optbench.df_to_arrays(prot)
ref = np.load(optbench.REF_PROT, allow_pickle=True)
rids, rvals = ref["ids"], ref["vals"]
assert np.array_equal(ids, rids)
d = np.abs(vals - rvals); d[np.isnan(d)] = 0
row_max = d.max(axis=1)
diff_rows = np.where(row_max > 1e-12)[0]
print("n proteins differing:", len(diff_rows))
# ion counts per protein in normed_df
pn = normed.index.get_level_values(0).to_numpy()
import collections
counts = collections.Counter(pn)
diff_prot_ids = [ids[i].split("\t")[0] for i in diff_rows]
ion_counts = sorted(counts[p] for p in diff_prot_ids)
print("ion-count distribution of differing proteins:", ion_counts[:5], "...", ion_counts[-5:] if len(ion_counts)>5 else "")
print("min ions among differing:", min(ion_counts), " max:", max(ion_counts))
print("all differing have >100 ions:", all(c>100 for c in ion_counts))
