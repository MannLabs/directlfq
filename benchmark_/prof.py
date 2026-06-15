import cProfile, pstats, io
import directlfq.config as config, directlfq.utils as lfqutils
import directlfq.normalization as lfqnorm, directlfq.protein_intensity_estimation as lfqprot
import optbench

normed = optbench.prep()
spec = list(
    lfqprot.get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan(
        normed, 10, 1
    )
)[:2500]
pr = cProfile.Profile()
pr.enable()
lfqprot.get_list_with_sequential_processing(iter(spec))
pr.disable()
s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats("tottime")
ps.print_stats(18)
print(s.getvalue())
