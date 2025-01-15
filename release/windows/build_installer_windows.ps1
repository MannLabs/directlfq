# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./*.egg-info
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build_pyinstaller
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist_pyinstaller

# Creating the wheel
python -m build
# Make sure you include the required extra packages and always use the stable options!

# substitute X.Y.Z-devN with X.Y.Z.devN
$WHL_NAME = "directlfq-0.3.1-dev0-py3-none-any.whl" -replace '(\d+\.\d+\.\d+)-dev(\d+)', '$1.dev$2'
pip install "dist/$WHL_NAME[stable,gui-stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/directlfq.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
# pip install jinja2==3.0
