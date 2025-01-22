#!/bin/bash
conda activate directlfq
cd ratio_tests
pip install pytest
pip install nbmake==1.5.3
echo "Running ratio tests"
python -m pytest --nbmake systemtest_ratio_mixed_species1.ipynb --nbmake-timeout=1800
python -m pytest --nbmake systemtest_ratio_mixed_species2.ipynb --nbmake-timeout=1800
cd ..
conda deactivate