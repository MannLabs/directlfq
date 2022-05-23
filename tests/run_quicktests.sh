conda activate directlfq
python download_testfiles.py quicktest
jupyter nbconvert --execute --to ./sensitivity_tests/lfq_set.ipynb
conda deactivate
