#!/bin/bash
set -e -u

# Build the installer for MacOS.
# This script must be run from the root of the repository.

rm -rf dist build *.egg-info
rm -rf dist_pyinstaller build_pyinstaller

# Creating the wheel
python setup.py sdist bdist_wheel
pip install "dist/directlfq-0.3.1-dev0-py3-none-any.whl[stable,gui]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/directlfq.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y

