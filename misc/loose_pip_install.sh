conda create -n directlfq python=3.8 -y
conda activate directlfq
pip install -e '../.[development]'
directlfq
conda deactivate
