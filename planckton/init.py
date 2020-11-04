import os

from foyer.forcefield import Forcefield
import mbuild as mb
import numpy as np
import parmed as pmd

from planckton.utils import base_units


class Compound(mb.Compound):
    """ Wrapper class for mb.Compound"""

    def __init__(self, path_to_mol2):
        super(Compound, self).__init__()
        mb.load(path_to_mol2, compound=self)
        # Calculate mass of compound
        self.mass = np.sum([atom.mass for atom in self.to_parmed().atoms])

        # This helps to_parmed use residues to apply ff more quickly
        self.name = os.path.basename(path_to_mol2).split(".")[0]

        # We need to rename the atom types
        compound_pmd = pmd.load_file(path_to_mol2)
        for atom_pmd, atom_mb in zip(compound_pmd, self):
            atom_mb.name = "_{}".format(atom_pmd.type)


class Pack:
    def __init__(
        self,
        compound,
        n_compounds,
        density,
        ff_file="compounds/gaff.4fxml",
        remove_hydrogen_atoms=False,
    ):
        if not isinstance(compound, (list, set)):
            self.compound = [compound]
        else:
            self.compound = compound
        if n_compounds is not None and not isinstance(n_compounds, (list, set)):
            self.n_compounds = [n_compounds]
        else:
            self.n_compounds = n_compounds

        self.residues = [comp.name for comp in self.compound]
        self.density = density
        self.ff = Forcefield(ff_file)
        self.remove_hydrogen_atoms = remove_hydrogen_atoms
        self.L = self._calculate_L()

    def _remove_hydrogen(self):
        for subcompound in self.compound:
            for atom in subcompound.particles():
                if atom.name in ["_hc", "_ha", "_h1", "_h4"]:
                    # NOTE: May not be a comprehensive list of
                    # all hydrogen types.
                    subcompound.remove(atom)

    def pack(self, box_expand_factor=5):
        """
        Parameters
        ----------
        box_expand_factor : float
            Factor by which to expand the box before packing for faster
            packing. (default 5)

        Returns
        -------
        typed_system : ParmEd structure
            ParmEd structure of filled box
        """
        units = base_units.base_units()

        if self.remove_hydrogen_atoms:
            self._remove_hydrogen()

        L = (self.L * box_expand_factor)
        # Extra factor to make packing faster, will shrink it out
        box = mb.Box([L,L,L])
        system = mb.packing.fill_box(
            self.compound,
            n_compounds=self.n_compounds,
            box=box,
            overlap=0.2,
        )
        system.box = box
        pmd_system = system.to_parmed(residues=[self.residues])
        typed_system = self.ff.apply(
                pmd_system,
                assert_angle_params=False,
                assert_dihedral_params=False
                )
        return typed_system

    def _calculate_L(self):
        total_mass = np.sum(
            [n * c.mass for c, n in zip(self.compound, self.n_compounds)]
        )
        # Conversion from (amu/(g/cm^3))**(1/3) to ang
        L = (total_mass / self.density) ** (1 / 3) * 1.1841763
        L /= 10  # convert ang to nm
        return L

