# TODO remove with old release workflow
conda create -n directlfq_pip_test python=3.8 -y
conda activate directlfq_pip_test
pip install "directlfq[stable, gui]"
directlfq
conda deactivate
