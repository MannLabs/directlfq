conda create -n directlfq python=3.9 -y
conda activate directlfq
pip install -e '../.[stable,development-stable, gui-stable]'
directlfq
conda deactivate
