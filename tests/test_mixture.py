from os import path, remove

from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


def test_mixture():
    pcbm = Compound(COMPOUND_FILE["PCBM"])
    p3ht = Compound(COMPOUND_FILE["P3HT"])
    packer = Pack(
        [pcbm, p3ht],
        ff_file=FORCE_FIELD["opv_gaff"],
        n_compounds=[2,3],
        density=0.01,
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
            shrink_time=1e3,
            )
    my_sim.run()


if __name__ == "__main__":
    if path.isfile("restart.gsd"):
        remove("restart.gsd")
    test_hydrogen_removal()
