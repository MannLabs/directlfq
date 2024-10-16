conda create -n directlfq python=3.9 -y
conda activate directlfq
pip install -e '../.[development-stable, gui]'
directlfq
conda deactivate
