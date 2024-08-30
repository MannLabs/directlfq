#!/bin/bash
set -u -e

# Build the install package for MacOS.
# This script must be run from the root of the repository after running build_installer_macos.sh

PACKAGE_NAME=directlfq
# BUILD_NAME is taken from environment variables, e.g. peptdeep-1.2.1-macos-darwin-arm64 or peptdeep-1.2.1-macos-darwin-x64
rm -rf ${BUILD_NAME}.pkg

# If needed, include additional source such as e.g.:
# cp ../../directlfq/data/*.fasta dist/directlfq/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/${PACKAGE_NAME}/Contents/Resources
cp release/logos/alpha_logo.icns dist/${PACKAGE_NAME}/Contents/Resources
mv dist/directlfq_gui dist/${PACKAGE_NAME}/Contents/MacOS
cp Info.plist dist/${PACKAGE_NAME}/Contents
cp directlfq_terminal dist/${PACKAGE_NAME}/Contents/MacOS
cp ./LICENSE Resources/LICENSE
cp release/logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/${PACKAGE_NAME} --identifier de.mpg.biochem.${PACKAGE_NAME}.app --version 0.2.20 --install-location /Applications/${PACKAGE_NAME}.app --scripts scripts ${PACKAGE_NAME}.pkg
productbuild --distribution distribution.xml --resources Resources --package-path ${PACKAGE_NAME}.pkg dist/${BUILD_NAME}.pkg
