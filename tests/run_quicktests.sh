conda activate directlfq
pip install pytest
pip install nbmake==1.5.3
cd quicktests
pytest --nbmake run_pipeline_w_different_input_formats.ipynb
directlfq lfq -i ../../test_data/system_tests/quicktests/diann/shortened_input.tsv
cd ..
conda deactivate