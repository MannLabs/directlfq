import numpy as np, optbench
import directlfq.protein_intensity_estimation as lfqprot
import pandas as pd
normed = optbench.prep()
pn = normed.index.get_level_values(0).to_numpy()
# find a protein with ~102 ions
import collections
counts = collections.Counter(pn)
target = [p for p,c in counts.items() if c==102][0]
sub = normed[pn==target]
vals = sub.to_numpy()
ion_names = sub.index.get_level_values(1).to_numpy()
# ProtvalCutter order
pc = lfqprot.ProtvalCutter(sub.copy(), maximum_df_length=100)
pc_idx = pc._sorted_idx[:100]
pc_ions = [t[1] for t in pc_idx]
# my lexsort order
with np.errstate(all="ignore"):
    neg = -np.nansum(vals, axis=1)
nanc = np.isnan(vals).sum(axis=1)
order = np.lexsort((neg, nanc))[:100]
my_ions = list(ion_names[order])
print("same set:", set(pc_ions)==set(my_ions))
print("same order:", pc_ions==my_ions)
# first divergence
for i,(a,b) in enumerate(zip(pc_ions, my_ions)):
    if a!=b:
        print(f"first divergence at {i}: pc={a} mine={b}")
        # show their keys
        ia = list(ion_names).index(a); ib = list(ion_names).index(b)
        print(f"  pc row {a}: nan={nanc[ia]} negsum={neg[ia]!r}")
        print(f"  my row {b}: nan={nanc[ib]} negsum={neg[ib]!r}")
        break
