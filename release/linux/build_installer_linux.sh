#!/bin/bash

# Build the installer for Linux.
# This script must be run from the root of the repository.

rm -rf dist
rm -rf build

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/directlfq-0.2.20-py3-none-any.whl[stable, gui]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/directlfq.spec -y
