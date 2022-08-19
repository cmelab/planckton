from os import path, remove

import unyt as u

from planckton.compounds import COMPOUND
from planckton.forcefields import FORCEFIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


def test_mixture():
    pcbm = Compound(COMPOUND["PCBM-gaff"])
    p3ht = Compound(COMPOUND["P3HT-gaff"])
    packer = Pack(
        [pcbm, p3ht],
        ff=FORCEFIELD["gaff-custom"],
        n_compounds=[2, 3],
        density=0.01 * u.g / u.cm**3,
    )

    system = packer.pack()
    my_sim = Simulation(
        system,
        kT=[3.0],
        tau=[1.0],
        n_steps=[1e3],
        gsd_write=1e2,
        log_write=1e2,
        shrink_steps=1e3,
    )
    my_sim.run()


if __name__ == "__main__":
    if path.isfile("restart.gsd"):
        remove("restart.gsd")
    test_mixture()
