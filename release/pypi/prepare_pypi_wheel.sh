# TODO remove with old release workflow
cd ../..
conda create -n directlfq_pypi_wheel python=3.8
conda activate directlfq_pypi_wheel
pip install twine
rm -rf dist
rm -rf build
python setup.py sdist bdist_wheel
twine check dist/*
conda deactivate
