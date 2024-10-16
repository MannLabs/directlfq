#!/bin/bash
conda activate directlfq
cd ratio_tests
pip install pytest
pip install nbmake==1.5.3
echo "Running ratio tests"
cd tests/ratio_tests
pytest --nbmake systemtest_ratio_mixed_species1.ipynb
pytest --nbmake systemtest_ratio_mixed_species2.ipynb
cd ../..
conda deactivate