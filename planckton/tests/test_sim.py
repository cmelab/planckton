import pytest
import unyt as u
from base_test import BaseTest

from planckton.compounds import COMPOUND
from planckton.forcefields import FORCEFIELD
from planckton.init import Compound, Pack
from planckton.sim import Simulation


class TestSimulations(BaseTest):
    @pytest.mark.parametrize("compound_name", COMPOUND.keys())
    def test_simple_sim(self, compound_name):
        compound = Compound(COMPOUND[compound_name])
        packer = Pack(compound, n_compounds=2, density=0.01 * u.g / u.cm**3)
        system = packer.pack()
        my_sim = Simulation(
            system,
            kT=[3.0],
            tau=[1.0],
            n_steps=[1e3],
            e_factor=0.5,
            gsd_write=1e2,
            log_write=1e2,
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
        )
        my_sim.run()

    def test_stack(self):
        compound = Compound(COMPOUND["P3HT-16-gaff"])
        packer = Pack(compound, n_compounds=2, density=0.01 * u.g / u.cm**3)
        system = packer.stack(separation=1.0)
        my_sim = Simulation(
            system,
            kT=[3.0],
            tau=[1.0],
            n_steps=[1e3],
            e_factor=0.5,
            gsd_write=1e2,
            log_write=1e2,
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
        )
        my_sim.run()

    def test_smiles_gaff(self):
        p3ht = Compound("c1cscc1CCCCCC")
        packer = Pack(
            p3ht,
            ff=FORCEFIELD["gaff"],
            n_compounds=2,
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
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
        )

    def test_gaff_noH(self):
        p3ht = Compound("c1cscc1CCCCCC")
        packer = Pack(
            p3ht,
            ff=FORCEFIELD["gaff"],
            n_compounds=2,
            density=0.01 * u.g / u.cm**3,
            remove_hydrogen_atoms=True,
        )
        system = packer.pack()
        my_sim = Simulation(
            system,
            kT=[3.0],
            tau=[1.0],
            n_steps=[1e3],
            gsd_write=1e2,
            log_write=1e2,
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
        )

    def test_smiles_opvgaff_raises(self):
        p3ht = Compound("c1cscc1CCCCCC")
        with pytest.raises(NotImplementedError):
            packer = Pack(
                p3ht,
                ff=FORCEFIELD["gaff-custom"],
                n_compounds=2,
                density=0.01 * u.g / u.cm**3,
            )

    def test_typed_gaff_raises(self):
        p3ht = Compound(COMPOUND["P3HT-16-gaff"])
        with pytest.warns(UserWarning):
            packer = Pack(
                p3ht,
                ff=FORCEFIELD["gaff"],
                n_compounds=2,
                density=0.01 * u.g / u.cm**3,
            )

    def test_nlist(self):
        p3ht = Compound("c1cscc1CCCCCC")
        packer = Pack(
            p3ht,
            ff=FORCEFIELD["gaff"],
            n_compounds=2,
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
            mode="cpu",
            shrink_steps=1e3,
            target_length=packer.L,
            nlist="tree",
        )

    def test_temps_ramp(self):
        p3ht = Compound("c1cscc1CCCCCC")
        packer = Pack(
            p3ht,
            ff=FORCEFIELD["gaff"],
            n_compounds=2,
            density=0.01 * u.g / u.cm**3,
        )
        system = packer.pack()
        my_sim = Simulation(
            system,
            kT=[3.0, 4.0],
            tau=[1.0, 1.0],
            n_steps=[1e3, 1e3],
            shrink_steps=1e3,
            target_length=packer.L,
            mode="cpu",
        )

    def test_bad_temps_raises(self):
        p3ht = Compound("c1cscc1CCCCCC")
        packer = Pack(
            p3ht,
            ff=FORCEFIELD["gaff"],
            n_compounds=2,
            density=0.01 * u.g / u.cm**3,
        )
        system = packer.pack()
        with pytest.raises(AssertionError):
            my_sim = Simulation(
                system,
                kT=[3.0, 4.0],
                tau=[1.0],
                n_steps=[1e3],
                mode="cpu",
            )
