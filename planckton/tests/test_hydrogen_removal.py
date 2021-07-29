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

    system = packer.pack()
    assert 1 not in [a.atomic_number for a in system.atoms]


def test_hydrogen_removal_and_sim():
    pcbm = Compound(COMPOUND["PCBM-gaff"])
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


def test_hydrogen_remove_gaff():
    p3ht = Compound("c1cscc1CCCCCC")
    p3ht_Hs = [h for h in p3ht.particles_by_element("H")]
    packer = Pack(
        p3ht,
        ff=FORCEFIELD["gaff"],
        n_compounds=2,
        density=0.01 * u.g / u.cm ** 3,
        remove_hydrogen_atoms=True,
    )
    system = packer.pack()
    assert "H" not in [a.name for a in system.atoms]
    assert p3ht.n_particles * 2 - len(system.atoms) == len(p3ht_Hs) * 2


if __name__ == "__main__":
    if path.isfile("restart.gsd"):
        remove("restart.gsd")
    test_hydrogen_removal()
    test_hydrogen_removal_and_sim()
