import pytest

from base_test import BaseTest
from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


class TestSimulations(BaseTest):
    @pytest.mark.parametrize("compound_name", COMPOUND_FILE.keys())
    def test_simple_sim(self, compound_name):
        compound = Compound(COMPOUND_FILE[compound_name])
        packer = Pack(
            compound,
            ff_file=FORCE_FIELD["opv_gaff"],
            n_compounds=2,
            density=0.01
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
            target_length=packer.L*10 # nm to A conversion
        )
        my_sim.run()
