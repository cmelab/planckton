from os import path, remove

import unyt as u

from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


def test_hydrogen_removal():
    pcbm = Compound(COMPOUND_FILE["PCBM"])
    packer = Pack(
        pcbm,
        ff=FORCE_FIELD["opv_gaff"],
        n_compounds=2,
        density=0.1 * u.g / u.cm**3,
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
    pcbm = Compound(COMPOUND_FILE["PCBM"])
    packer = Pack(
        pcbm,
        ff=FORCE_FIELD["opv_gaff"],
        n_compounds=2,
        density=0.1 * u.g / u.cm**3,
        remove_hydrogen_atoms=True,
    )
    print("packer init")
    system = packer.pack()
    print("packer packed")
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
    print("sim init")
    my_sim.run()


if __name__ == "__main__":
    if path.isfile("restart.gsd"):
        remove("restart.gsd")
    test_hydrogen_removal()
    test_hydrogen_removal_and_sim()
