conda activate directlfq
pip install wget
python download_testfiles.py quicktest
cd quicktests
jupyter nbconvert --to script run_pipeline_w_different_input_formats.ipynb
python run_pipeline_w_different_input_formats.py
cd ..
conda deactivate