conda activate directlfq
pip install wget
python download_testfiles.py quicktest
jupyter nbconvert --NotebookApp.kernel_manager_class=notebook.services.kernels.kernelmanager.AsyncMappingKernelManager --to notebook --execute ./quicktests/run_pipeline_w_different_input_formats.ipynb
conda deactivate
