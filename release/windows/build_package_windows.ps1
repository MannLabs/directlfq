# Build the install package for Windows.
# This script must be run from the root of the repository after running build_installer_windows.ps1


# Apparently, ISCC uses the directory of the .iss input file as the working directory.
mkdir release/windows/dist
mv dist_pyinstaller/* release/windows/dist

# Wrapping the pyinstaller folder in a .exe package
&  "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" .\release\windows\directlfq_innoinstaller.iss
# release workflow expects artifact at root of repository
mv .\release\windows\dist\*.exe .
