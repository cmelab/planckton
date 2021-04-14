import pytest

import unyt as u
from unyt.exceptions import UnitConversionError

from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from planckton.init import Compound, Pack


def test_load_smiles():
    p3ht = Compound("c1cscc1CCCCCC")


def test_bad_units():
    pcbm = Compound(COMPOUND_FILE["PCBM"])
    with pytest.raises(UnitConversionError):
        packer = Pack(
            pcbm, ff=FORCE_FIELD["opv_gaff"], n_compounds=2, density=2 * u.m
        )


def test_no_units():
    pcbm = Compound(COMPOUND_FILE["PCBM"])
    with pytest.raises(TypeError):
        packer = Pack(
            pcbm, ff=FORCE_FIELD["opv_gaff"], n_compounds=2, density=2
        )
