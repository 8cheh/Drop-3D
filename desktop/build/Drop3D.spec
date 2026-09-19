# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller specification for the Drop-3D desktop application.

Build with ``desktop/build/Build-Windows.ps1``, which creates a clean virtual
environment first.  Running PyInstaller from an environment that can see a shared
site-packages is the one mistake that matters here: it pulls in every package the
machine happens to have, and the first attempt at this produced a 932 MB bundle
containing torch, OpenCV and transformers -- none of which this application
imports.  Frozen from a clean environment it is about 126 MB.

The interface files are shipped as data.  Without them the executable starts,
fails to find index.html and exits into a native error dialog, which is a
confusing way to learn that ``--add-data`` was forgotten.
"""

import os

from PyInstaller.utils.hooks import collect_all

SPEC_DIR = os.path.abspath(SPECPATH)                    # desktop/build
DESKTOP = os.path.dirname(SPEC_DIR)                     # desktop
SRC = os.path.join(DESKTOP, 'src')
ASSETS = os.path.join(SRC, 'drop3d_desktop', 'assets')

datas = [(ASSETS, 'assets')]
binaries = []
hiddenimports = []

# pywebview and its .NET bridge load parts of themselves at runtime rather than
# through an import statement, so their whole package trees have to be collected
# explicitly or the frozen application fails on first use.
for package in ('webview', 'clr_loader', 'pythonnet'):
    pkg_datas, pkg_binaries, pkg_hidden = collect_all(package)
    datas += pkg_datas
    binaries += pkg_binaries
    hiddenimports += pkg_hidden

hiddenimports += [
    'drop3d',
    'drop3d_desktop',
    'drop3d.tools.shape_parameter_scan',
]

# Nothing here is imported by the application.  matplotlib is installed on many
# development machines and drags in a large tree; the GUI toolkits are pywebview
# alternatives that are not used.
excludes = [
    'matplotlib', 'tkinter', 'PyQt5', 'PyQt6', 'PySide2', 'PySide6',
    'IPython', 'jupyter', 'pytest', 'setuptools._distutils',
]

a = Analysis(
    [os.path.join(SPEC_DIR, 'launch.py')],
    pathex=[SRC],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='Drop-3D',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,          # UPX on scipy's DLLs is a known source of corruption
    console=False,      # no console window; startup errors go to the log file
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name='Drop-3D',
)
