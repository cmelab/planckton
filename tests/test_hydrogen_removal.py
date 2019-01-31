from planckton.init import Compound, Pack
from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from os import path, remove


def test_hydrogen_removal():
    pcbm = Compound(COMPOUND_FILE["PCBM"])
    packer = Pack(
        pcbm,
        ff_file=FORCE_FIELD["opv_gaff"],
        n_compounds=2,
        density=0.01,
        out_file="test_init.hoomdxml",
        remove_hydrogen_atoms=True,
    )

    packer._remove_hydrogen()

    for atom in packer.compound[0].particles():
        assert atom.name not in [
            "_hc",
            "_ha",
            "_h1",
            "_h4",
        ], "Hydrogen found in system!"


if __name__ == "__main__":
    if path.isfile("restart.gsd"):
        remove("restart.gsd")
    test_hydrogen_removal()
