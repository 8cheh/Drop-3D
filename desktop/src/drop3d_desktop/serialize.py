"""JSON-safe conversion for the JavaScript bridge.

Everything that crosses into the web view has to be plain JSON.  ``drop3d``'s
public API returns dataclasses and dictionaries, which are *almost* JSON already,
so this module is the single place that knows how to flatten them.  Keeping that
knowledge in one file means a new result type only has to be taught once, and the
bridge functions stay readable.

Two things are deliberately not JSON and need handling here:

* ``numpy`` scalars and arrays -- ``float64`` is not ``float`` and ``json`` will
  refuse it;
* ``NaN`` / ``Infinity`` -- the JSON spec has no such literals.  ``json.dumps``
  would happily emit them and the browser's ``JSON.parse`` would then throw,
  which is exactly the kind of failure that is hard to trace back to its cause.
  They are reported as ``null`` instead, so the interface can show "not
  available" rather than dying.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import math
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ['to_jsonable', 'dumps', 'round_floats']


def _clean_scalar(value: Any) -> Any:
    """Make one non-container value JSON-safe."""
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    elif isinstance(value, np.ndarray):
        return [_clean_scalar(v) for v in value.tolist()]
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, (Path, dt.date, dt.datetime)):
        return str(value)
    return value


def _dataclass_to_dict(obj: Any, depth: int, include_properties: bool) -> dict[str, Any]:
    """Flatten a dataclass, including the computed properties it exposes.

    Skipping the properties would be the obvious implementation and it is the
    wrong one.  ``drop3d`` puts its most load-bearing summaries in properties
    rather than fields -- ``ValidityReport.error_codes`` and
    ``reliability_class``, ``DynamicsReport.verdict``, ``EllipseFit.aspect`` --
    so a fields-only walk would hand the interface a gate report with no verdict
    and no error codes, and the silence would look like a missing feature rather
    than a serialisation bug.

    Properties are read on a best-effort basis: one that raises (or needs an
    argument) is skipped instead of taking the whole payload down.
    """
    data = {
        f.name: to_jsonable(getattr(obj, f.name), depth + 1, include_properties)
        for f in dataclasses.fields(obj)
    }
    if include_properties:
        for name, attr in vars(type(obj)).items():
            if name.startswith('_') or name in data or not isinstance(attr, property):
                continue
            try:
                data[name] = to_jsonable(getattr(obj, name), depth + 1, include_properties)
            except Exception:  # noqa: BLE001 - a property is a convenience, not a contract
                continue
    return data


def to_jsonable(obj: Any, _depth: int = 0, include_properties: bool = True) -> Any:
    """Recursively convert ``obj`` into something ``json.dumps`` accepts.

    The depth guard exists because a malformed or self-referential object would
    otherwise recurse until the interpreter's stack runs out -- inside a GUI that
    surfaces as a frozen window rather than a visible error.
    """
    if _depth > 12:
        return None

    if obj is None or isinstance(obj, (str, bool, int)):
        return obj
    if isinstance(obj, (float, np.floating, np.integer, np.bool_)):
        return _clean_scalar(obj)
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return _dataclass_to_dict(obj, _depth, include_properties)
    if isinstance(obj, dict):
        return {
            str(k): to_jsonable(v, _depth + 1, include_properties)
            for k, v in obj.items()
        }
    if isinstance(obj, (list, tuple, set, frozenset)):
        return [to_jsonable(v, _depth + 1, include_properties) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist(), _depth + 1, include_properties)
    if isinstance(obj, (Path, dt.date, dt.datetime)):
        return str(obj)
    # Anything else (a solver object, a callable, an image handle) is not
    # transportable.  Saying so is better than dropping the key silently, which
    # would make the interface show a blank where a number was expected.
    return f'<{type(obj).__name__} not serialisable>'


def dumps(obj: Any) -> str:
    """``json.dumps`` that will not emit NaN/Infinity.

    ``allow_nan=False`` is the point: if a NaN slips past :func:`to_jsonable`,
    this raises here, in Python, where the traceback names the field -- rather
    than in the browser as an unhelpful parse error.
    """
    import json

    return json.dumps(to_jsonable(obj), ensure_ascii=False, allow_nan=False)


def round_floats(obj: Any, digits: int = 6) -> Any:
    """Round every float in a JSON-able tree.

    The bridge payloads carry residual arrays and scanned curves; full float64
    noise in those would bloat the message for no visible gain.
    """
    if isinstance(obj, float):
        return round(obj, digits)
    if isinstance(obj, dict):
        return {k: round_floats(v, digits) for k, v in obj.items()}
    if isinstance(obj, list):
        return [round_floats(v, digits) for v in obj]
    return obj
