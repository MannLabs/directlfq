#!/bin/bash
conda activate directlfq
cd ratio_tests
pip install pytest
pip install nbmake==1.5.3
echo "Running ratio tests"
pytest --nbmake systemtest_ratio_mixed_species1.ipynb --nbmake-timeout=1200
pytest --nbmake systemtest_ratio_mixed_species2.ipynb --nbmake-timeout=1200
cd ..
conda deactivate