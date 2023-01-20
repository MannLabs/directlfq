conda activate directlfq
pip install wget
python download_testfiles.py quicktest
cd quicktests
echo "Script executed from: ${PWD}"
papermill run_pipeline_w_different_input_formats.ipynb run_pipeline_w_different_input_formats.nbconvert.ipynb
cd ..
conda deactivate