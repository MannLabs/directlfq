import numpy as np, optbench, collections
import directlfq.protein_intensity_estimation as lfqprot

normed = optbench.prep()
pn = normed.index.get_level_values(0).to_numpy()
# differing proteins from before: recompute
prot, ion = lfqprot.estimate_protein_intensities(normed, 1, 10, 1)
ids, vals_out = optbench.df_to_arrays(prot)
ref = np.load(optbench.REF_PROT, allow_pickle=True)
d = np.abs(vals_out - ref["vals"])
d[np.isnan(d)] = 0
diff_ids = [ids[i].split("\t")[0] for i in np.where(d.max(axis=1) > 1e-12)[0]]
print("checking", len(diff_ids), "differing proteins")
for p in diff_ids[:4]:
    sub = normed[pn == p]
    # OLD path: ProtvalCutter df -> nansum(2**df)
    cut_df = optbench.protvalcutter_reference_cut(sub, 100)
    old_sum = np.nansum(2**cut_df)
    # NEW path: array cut -> nansum(2**arr)
    vals = sub.to_numpy()
    ino = sub.index.get_level_values(1).to_numpy()
    nv, ni = lfqprot._cut_peptide_values(vals, ino, 100)
    new_sum = np.nansum(2**nv)
    # also compare order
    order_same = [t[1] for t in cut_df.index] == list(ni)
    print(
        f"{p[:20]:22s} old_sum={old_sum!r} new_sum={new_sum!r} equal={old_sum == new_sum} order_same={order_same}"
    )
