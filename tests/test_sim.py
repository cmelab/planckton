from planckton.init import Compound, Pack
from planckton.sim import Simulation
from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD


def test_simple_sim():
    pcbm = Compound(COMPOUND_FILE["pcbm"])
    packer = Pack(pcbm, ff_file=FORCE_FIELD["opv_gaff"], n_compounds=2, density=0.01)
    packer.pack()
    my_sim = Simulation(
        "init.hoomdxml", kT=3.0, gsd_write=1e2, log_write=1e2, e_factor=0.5, n_steps=3e3, mode="gpu", shrink_time=1e3
    )
    my_sim.run()

if __name__ == "__main__":
    test_simple_sim()
