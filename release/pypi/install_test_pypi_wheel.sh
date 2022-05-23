conda create -n directlfq_pip_test python=3.8 -y
conda activate directlfq_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "directlfq[stable]"
directlfq
conda deactivate
