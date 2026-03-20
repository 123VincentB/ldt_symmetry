# -*- coding: utf-8 -*-
"""ldt_symmetry
===============

Extension to `pyldt <https://pypi.org/project/eulumdat-py/>`_ for
symmetrising EULUMDAT (.ldt) photometric files.

Public API
----------
:class:`LdtSymmetriser`
    Symmetrise a ``pyldt.Ldt`` object for a given ISYM mode (0–4).

:class:`LdtAutoDetector`
    Automatically detect the most appropriate ISYM mode.
    Optional — use only when the symmetry mode is not known in advance.

Quick start
-----------
::

    from pyldt import LdtReader, LdtWriter
    from ldt_symmetry import LdtSymmetriser

    ldt = LdtReader.read("luminaire.ldt")
    ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2)
    LdtWriter.write(ldt_sym, "luminaire_SYM.ldt")

With automatic detection::

    from pyldt import LdtReader, LdtWriter
    from ldt_symmetry import LdtAutoDetector, LdtSymmetriser

    ldt = LdtReader.read("luminaire.ldt")
    detector = LdtAutoDetector()
    isym = detector.detect(ldt)
    ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
    LdtWriter.write(ldt_sym, "luminaire_SYM.ldt")
"""

from .auto_detector import LdtAutoDetector
from .symmetriser import LdtSymmetriser

__all__ = ["LdtSymmetriser", "LdtAutoDetector"]
__version__ = "1.0.0"
