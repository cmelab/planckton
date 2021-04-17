from os import path, remove

import unyt as u

from planckton.compounds import COMPOUND
from planckton.forcefields import FORCEFIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


def test_hydrogen_removal():
    pcbm = Compound(COMPOUND["PCBM-gaff"])
    packer = Pack(
        pcbm,
        ff=FORCEFIELD["gaff-custom"],
        n_compounds=2,
        density=0.1 * u.g / u.cm ** 3,
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


def test_hydrogen_removal_and_sim():
    pcbm = Compound(COMPOUND_FILE["PCBM-gaff"])
    packer = Pack(
        pcbm,
        ff=FORCEFIELD["gaff-custom"],
        n_compounds=2,
        density=0.1 * u.g / u.cm ** 3,
        remove_hydrogen_atoms=True,
    )
    system = packer.pack()
    my_sim = Simulation(
        system,
        kT=3.0,
        gsd_write=1e2,
        log_write=1e2,
        e_factor=0.5,
        n_steps=3e3,
        mode="cpu",
        shrink_steps=1e3,
    )
    my_sim.run()


if __name__ == "__main__":
    if path.isfile("restart.gsd"):
        remove("restart.gsd")
    test_hydrogen_removal()
    test_hydrogen_removal_and_sim()
