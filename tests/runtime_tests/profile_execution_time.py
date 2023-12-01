input_file = "/Users/constantin/workspace/directlfq/test_data/system_tests/Spectronaut_LargeFC/re_run_SN15/20220607_153923_MP-LFC-MS1var-OT-S1-120kMS1_Report.tsv.spectronaut_fragion_isotopes.aq_reformat.tsv"


import directlfq.lfq_manager as lfqmgr

lfqmgr.run_lfq(input_file=input_file, num_cores=1)