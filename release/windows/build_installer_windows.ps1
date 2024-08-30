# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist

# Creating the wheel
python setup.py sdist bdist_wheel
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/directlfq-0.2.20-py3-none-any.whl[stable, gui]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/directlfq.spec -y
pip install jinja2==3.0

# for some reason, the installer builder expects the files here
mkdir release/windows/dist
mv dist/* release/windows/dist
mkdir release/windows/dist/directlfq_gui
mv release/windows/dist/directlfq_gui.exe release/windows/dist/directlfq_gui

