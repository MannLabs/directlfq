conda activate directlfq
pip install pytest
pip install nbmake==1.5.3
pip install alphabase>=1.2.1
echo "Current working directory before cd $(pwd)"
cd quicktests
echo "Current working directory after cd $(pwd)"
pytest --nbmake run_pipeline_w_different_input_formats.ipynb
echo "Current working directory after cd and nb execution $(pwd)"
directlfq lfq -i ../../test_data/system_tests/quicktests/diann/shortened_input.tsv
cd ..
conda deactivate