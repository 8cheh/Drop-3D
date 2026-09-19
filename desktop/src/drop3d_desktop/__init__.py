"""Desktop interface for Drop-3D.

The application is two halves that meet at one bridge:

* this Python side, which owns the window and calls ``drop3d``;
* the web view, which owns everything the user sees.

They communicate only through :class:`drop3d_desktop.api.DesktopApi`.  No HTTP
server, no localhost port, no open socket -- the bridge is an in-process object,
which means the packaged application has no network surface at all.
"""

from __future__ import annotations

__version__ = '0.1.0'

__all__ = ['__version__']
