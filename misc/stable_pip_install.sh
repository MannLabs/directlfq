conda create -n directlfq python=3.9 -y
conda activate directlfq
pip install -e '../.[stable, gui-stable, development]'
directlfq
conda deactivate
