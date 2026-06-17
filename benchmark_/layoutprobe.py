import numpy as np, optbench, collections
import directlfq.protein_intensity_estimation as lfqprot

normed = optbench.prep()
pn = normed.index.get_level_values(0).to_numpy()
prot, _ = lfqprot.estimate_protein_intensities(normed, 1, 10, 1)
ids, vo = optbench.df_to_arrays(prot)
d = np.abs(vo - np.load(optbench.REF_PROT, allow_pickle=True)["vals"])
d[np.isnan(d)] = 0
diff_ids = [ids[i].split("\t")[0] for i in np.where(d.max(axis=1) > 1e-12)[0]]
ok = {"C": 0, "F": 0, "T": 0, "col": 0}
for p in diff_ids:
    sub = normed[pn == p]
    cut_df = optbench.protvalcutter_reference_cut(sub.copy(), 100)
    old = np.nansum(2**cut_df)
    vals = sub.to_numpy()
    ino = sub.index.get_level_values(1).to_numpy()
    nv, _ = lfqprot._cut_peptide_values(vals, ino, 100)
    pw = 2**nv
    if np.nansum(pw) == old:
        ok["C"] += 1
    if np.nansum(np.asfortranarray(pw)) == old:
        ok["F"] += 1
    if np.nansum(pw.T) == old:
        ok["T"] += 1
    if np.nansum(np.nansum(pw, axis=0)) == old:
        ok["col"] += 1
print("matches old_sum (out of %d):" % len(diff_ids), ok)
