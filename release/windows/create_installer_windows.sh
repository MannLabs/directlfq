#!bash

# TODO remove with old release workflow

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n directlfq_installer python=3.9 -y
conda activate directlfq_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/windows
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/directlfq-0.3.0-py3-none-any.whl[stable, gui]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/directlfq.spec -y
pip install jinja2==3.0
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../directlfq/data/*.fasta dist/directlfq/data

# Wrapping the pyinstaller folder in a .exe package
"C:\Program Files (x86)\Inno Setup 6\ISCC.exe" directlfq_innoinstaller_old.iss
# WARNING: this assumes a static location for innosetup
