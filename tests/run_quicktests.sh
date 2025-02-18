conda activate directlfq
pip install -r ../requirements/requirements_development.txt
cd quicktests
python -m pytest --nbmake run_pipeline_w_different_input_formats.ipynb
directlfq lfq -i ../../test_data/system_tests/quicktests/diann/shortened_input.tsv
cd ..
conda deactivate