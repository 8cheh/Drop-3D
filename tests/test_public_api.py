"""Guards on the package's public surface.

``__all__`` is a promise: it is what ``from drop3d import *`` resolves, what the
documentation indexes, and what a reader trusts when they look up a name.  A name
listed there but never imported breaks the star import outright and is invisible
to every other test in this suite, because no test used that name.

This was not hypothetical.  ``footprint_from_contact_line`` was listed in
``__all__`` and never imported, so ``from drop3d import *`` raised

    AttributeError: module 'drop3d' has no attribute 'footprint_from_contact_line'

while the function itself worked fine when imported from its own module -- which
is exactly the shape of bug that survives a green test suite.
"""

from __future__ import annotations

import drop3d


def test_every_advertised_name_exists():
    missing = sorted(name for name in drop3d.__all__ if not hasattr(drop3d, name))
    assert not missing, (
        f'{missing} are listed in drop3d.__all__ but are not importable from the '
        f'package; add them to the import block in src/drop3d/__init__.py'
    )


def test_star_import_works():
    """The consequence of the bug above, tested directly."""
    namespace: dict = {}
    exec('from drop3d import *', namespace)  # noqa: S102 - the point is the import
    assert 'footprint_from_contact_line' in namespace


def test_all_has_no_duplicates():
    duplicates = sorted({name for name in drop3d.__all__ if drop3d.__all__.count(name) > 1})
    assert not duplicates, f'duplicated in __all__: {duplicates}'


def test_every_advertised_name_is_a_callable_or_a_constant():
    """Catches a name that exists only because a submodule happened to be
    imported, which would make the public surface depend on import order."""
    for name in drop3d.__all__:
        value = getattr(drop3d, name)
        assert not hasattr(value, '__path__'), (
            f'{name} is a module, not a public symbol; import it as drop3d.{name}'
        )
