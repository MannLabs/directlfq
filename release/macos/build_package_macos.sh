#!/bin/bash
set -u -e

# Build the install package for MacOS.
# This script must be run from the root of the repository after running build_installer_macos.sh

PACKAGE_NAME=directlfq
# BUILD_NAME is taken from environment variables, e.g. directlfq-1.2.3-macos-darwin-arm64 or directlfq-1.2.3-macos-darwin-x64
rm -rf ${BUILD_NAME}.pkg

# If needed, include additional source such as e.g.:
# cp ../../directlfq/data/*.fasta dist/directlfq/data

# Wrapping the pyinstaller folder in a .pkg package
CONTENTS_FOLDER=dist_pyinstaller/${PACKAGE_NAME}/Contents

mkdir -p ${CONTENTS_FOLDER}/Resources
cp release/logos/alpha_logo.icns ${CONTENTS_FOLDER}/Resources
mv dist_pyinstaller/directlfq_gui ${CONTENTS_FOLDER}/MacOS
cp release/macos/Info.plist ${CONTENTS_FOLDER}
cp release/macos/directlfq_terminal ${CONTENTS_FOLDER}/MacOS
cp ./LICENSE ${CONTENTS_FOLDER}/Resources/LICENSE
cp release/logos/alpha_logo.png ${CONTENTS_FOLDER}/Resources

# link _internal folder containing the python libraries to the Frameworks folder where they are expected
# to avoid e.g. "Failed to load Python shared library '/Applications/AlphaMap.app/Contents/Frameworks/libpython3.8.dylib'"
cd ${CONTENTS_FOLDER}
ln -s ./MacOS/_internal ./Frameworks
cd -

chmod 777 release/macos/scripts/*

pkgbuild --root dist_pyinstaller/${PACKAGE_NAME} --identifier de.mpg.biochem.${PACKAGE_NAME}.app --version 0.3.3-dev0 --install-location /Applications/${PACKAGE_NAME}.app --scripts release/macos/scripts ${PACKAGE_NAME}.pkg
productbuild --distribution release/macos/distribution.xml --resources release/macos/Resources --package-path ${PACKAGE_NAME}.pkg ${BUILD_NAME}.pkg
