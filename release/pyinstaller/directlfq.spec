# -*- mode: python ; coding: utf-8 -*-

import os
import sys
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT, BUNDLE, TOC
import PyInstaller.utils.hooks


##################### User definitions
exe_name = 'directlfq_gui'
script_name = 'directlfq_pyinstaller.py'
if sys.platform[:6] == "darwin":
	icon = '../logos/alpha_logo.icns'
else:
	icon = '../logos/alpha_logo.ico'
block_cipher = None
location = os.getcwd()
project = "directlfq"
bundle_name = "directlfq"
#####################


datas, binaries, hidden_imports = PyInstaller.utils.hooks.collect_all(
	project,
	include_py_files=True
)

# in case we need other alphaX packages:
# alpha_x = ['alphabase', ]
# for alpha_package in alpha_x:
# 	_datas, _binaries, _hidden_imports = PyInstaller.utils.hooks.collect_all(
# 		alpha_package,
# 		include_py_files=True
# 	)
# 	datas+=_datas
# 	binaries+=_binaries
# 	hidden_imports+=_hidden_imports

hidden_imports = [h for h in hidden_imports if "__pycache__" not in h]
hidden_imports = sorted(
		[h for h in hidden_imports if "tests" not in h.split(".")]
	)
datas = [d for d in datas if ("__pycache__" not in d[0]) and (d[1] not in [".", "Resources", "scripts"])]

a = Analysis(
	[script_name],
	pathex=[location],
	binaries=binaries,
	datas=datas,
	hiddenimports=hidden_imports,
	hookspath=[],
	runtime_hooks=[],
	excludes=[h for h in hidden_imports if "datashader" in h],
	win_no_prefer_redirects=False,
	win_private_assemblies=False,
	cipher=block_cipher,
	noarchive=False
)
pyz = PYZ(
	a.pure,
	a.zipped_data,
	cipher=block_cipher
)

if sys.platform[:5] == "linux":
	exe = EXE(
		pyz,
		a.scripts,
		a.binaries,
		a.zipfiles,
		a.datas,
		name=bundle_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		upx_exclude=[],
		icon=icon
	)
else:
	exe = EXE(
		pyz,
		a.scripts,
		# a.binaries,
		a.zipfiles,
		# a.datas,
		exclude_binaries=True,
		name=exe_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		icon=icon
	)
	coll = COLLECT(
		exe,
		a.binaries,
		# a.zipfiles,
		a.datas,
		strip=False,
		upx=True,
		upx_exclude=[],
		name=exe_name
	)
