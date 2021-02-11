import pytest

from base_test import BaseTest
from planckton.compounds import COMPOUND_FILE
from planckton.force_fields import FORCE_FIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


class TestRigid(BaseTest):
    @pytest.mark.parametrize("compound_name", COMPOUND_FILE.keys())
    def test_rigid_bodies(self, compound_name):
        compound = Compound(COMPOUND_FILE[compound_name], rigid=True)
        packer = Pack(
            compound, ff=FORCE_FIELD["opv_gaff"], n_compounds=2, density=0.01
        )
        system = packer.pack()
        my_sim = Simulation(
            system,
            kT=3.0,
            gsd_write=1e2,
            log_write=1e2,
            rigid_inds=packer.rigid_inds,
            rigid_typeids=packer.rigid_typeids,
            n_steps=3e3,
            mode="cpu",
            shrink_time=1e3,
        )
        my_sim.run()
