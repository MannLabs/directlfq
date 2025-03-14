#!/bin/bash
set -e -u

# Build the install package for Linux.
# This script must be run from the root of the repository after running build_installer_linux.sh

PACKAGE_NAME=directlfq

# BUILD_NAME is taken from environment variables, e.g. 'directlfq-1.2.3-linux-x64'
rm -rf ${BUILD_NAME}.deb

# If needed, include additional source such as e.g.:
# cp ../../directlfq/data/*.fasta dist/directlfq/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mkdir -p dist_pyinstaller/${BUILD_NAME}/usr/local/bin
mv dist_pyinstaller/${PACKAGE_NAME} dist_pyinstaller/${BUILD_NAME}/usr/local/bin/${PACKAGE_NAME}
mkdir dist_pyinstaller/${BUILD_NAME}/DEBIAN
cp release/linux/control dist_pyinstaller/${BUILD_NAME}/DEBIAN
dpkg-deb --build --root-owner-group dist_pyinstaller/${BUILD_NAME}

# release workflow expects artifact at root of repository
mv dist_pyinstaller/${BUILD_NAME}.deb .
