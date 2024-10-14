#!/bin/bash
conda activate directlfq
pip install pytest
pip install nbmake==1.5.3
echo "Running ratio tests"
pytest --nbmake tests/ratio_tests/systemtest_ratio_mixed_species1.ipynb
pytest --nbmake tests/ratio_tests/systemtest_ratio_mixed_species2.ipynb
conda deactivate