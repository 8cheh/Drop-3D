"""Freeze entry point.

PyInstaller has to start from a plain script, and a plain script has no package
context -- so starting it on ``drop3d_desktop/__main__.py`` directly would break
that module's relative imports.  This one-line shim imports the package properly
instead, which keeps the application importable both frozen and from source.
"""

from drop3d_desktop.__main__ import main

if __name__ == '__main__':
    raise SystemExit(main())
