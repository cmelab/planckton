import pytest
import unyt as u
from unyt.exceptions import UnitConversionError

from planckton.compounds import COMPOUND
from planckton.init import Compound, Pack


def test_load_smiles():
    p3ht = Compound("c1cscc1CCCCCC")


def test_bad_units():
    pcbm = Compound(COMPOUND["PCBM-gaff"])
    with pytest.raises(UnitConversionError):
        packer = Pack(pcbm, n_compounds=2, density=2 * u.m)


def test_no_units():
    pcbm = Compound(COMPOUND["PCBM-gaff"])
    with pytest.raises(TypeError):
        packer = Pack(pcbm, n_compounds=2, density=2)
