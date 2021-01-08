from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from planckton.init import Compound, Pack


def test_load_smiles():
    p3ht = Compound("c1cscc1CCCCCC")
