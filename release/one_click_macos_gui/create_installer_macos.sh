#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=directlfq.pkg
if test -f "$FILE"; then
  rm directlfq.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n directlfqinstaller python=3.8 -y
conda activate directlfqinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/directlfq-0.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/directlfq.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../directlfq/data/*.fasta dist/directlfq/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/directlfq/Contents/Resources
cp ../logos/alpha_logo.icns dist/directlfq/Contents/Resources
mv dist/directlfq_gui dist/directlfq/Contents/MacOS
cp Info.plist dist/directlfq/Contents
cp directlfq_terminal dist/directlfq/Contents/MacOS
cp ../../LICENSE Resources/LICENSE
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/directlfq --identifier de.mpg.biochem.directlfq.app --version 0.2.1 --install-location /Applications/directlfq.app --scripts scripts directlfq.pkg
productbuild --distribution distribution.xml --resources Resources --package-path directlfq.pkg dist/directlfq_gui_installer_macos.pkg
