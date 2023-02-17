#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n directlfq_installer python=3.8 -y
conda activate directlfq_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_linux_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/directlfq-0.2.4-py3-none-any.whl[stable, gui]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/directlfq.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../directlfq/data/*.fasta dist/directlfq/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mkdir -p dist/directlfq_gui_installer_linux/usr/local/bin
mv dist/directlfq dist/directlfq_gui_installer_linux/usr/local/bin/directlfq
mkdir dist/directlfq_gui_installer_linux/DEBIAN
cp control dist/directlfq_gui_installer_linux/DEBIAN
dpkg-deb --build --root-owner-group dist/directlfq_gui_installer_linux/
