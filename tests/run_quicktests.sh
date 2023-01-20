conda activate directlfq
pip install wget
pip install papermill
python download_testfiles.py quicktest
cd quicktests
echo "Script executed from: ${PWD}"
jupyter execute run_pipeline_w_different_input_formats.ipynb
cd ..
conda deactivate