import pytest

from base_test import BaseTest
import unyt as u

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
            ff=FORCE_FIELD["opv_gaff"],
            n_compounds=2,
            density=0.01 * u.g / u.cm**3
        )
        system = packer.pack()
        my_sim = Simulation(
            system,
            kT=3.0,
            gsd_write=1e2,
            log_write=1e2,
            e_factor=1,
            n_steps=3e3,
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
        )
        my_sim.run()

    def test_smiles_gaff(self):
        p3ht = Compound("c1cscc1CCCCCC")
        packer = Pack(
            p3ht,
            ff=FORCE_FIELD["gaff"],
            n_compounds=2,
            density=0.01 * u.g / u.cm**3
        )
        system = packer.pack()
        my_sim = Simulation(
            system,
            kT=3.0,
            gsd_write=1e2,
            log_write=1e2,
            e_factor=1,
            n_steps=3e3,
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
        )

    def test_gaff_noH_raises(self):
        p3ht = Compound("c1cscc1CCCCCC")
        with pytest.raises(NotImplementedError):
            packer = Pack(
                p3ht,
                ff=FORCE_FIELD["gaff"],
                n_compounds=2,
                density=0.01 * u.g / u.cm**3,
                remove_hydrogen_atoms=True
            )

    def test_smiles_opvgaff_raises(self):
        p3ht = Compound("c1cscc1CCCCCC")
        with pytest.raises(NotImplementedError):
            packer = Pack(
                p3ht,
                ff=FORCE_FIELD["opv_gaff"],
                n_compounds=2,
                density=0.01 * u.g / u.cm**3,
            )

    def test_typed_gaff_raises(self):
        p3ht = Compound(COMPOUND_FILE["P3HT_16"])
        with pytest.raises(NotImplementedError):
            packer = Pack(
                p3ht,
                ff=FORCE_FIELD["gaff"],
                n_compounds=2,
                density=0.01 * u.g / u.cm**3,
            )
