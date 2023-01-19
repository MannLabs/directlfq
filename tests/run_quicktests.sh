conda activate directlfq
pip install wget
python download_testfiles.py quicktest
jupyter nbconvert --to notebook --execute ./quicktests/run_pipeline_w_different_input_formats.ipynb
conda deactivate
