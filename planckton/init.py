import os

import foyer
import mbuild as mb
import numpy as np
import parmed as pmd

from planckton.force_fields import FORCE_FIELD
from planckton.utils import base_units
from planckton.utils.rigid import connect_rings


class Compound(mb.Compound):
    """ Wrapper class for mb.Compound"""

    def __init__(self, input_str, rigid=False):
        super(Compound, self).__init__()
        if os.path.exists(input_str):
            mb.load(input_str, compound=self)
        else:
            mb.load(input_str, smiles=True, compound=self)

        # Calculate mass of compound
        self.mass = np.sum([atom.mass for atom in self.to_parmed().atoms])

        # This helps to_parmed use residues to apply ff more quickly
        self.name = os.path.basename(input_str).split(".")[0]

        if rigid:
            # Find conjugated rings and determine how many rigid bodies
            # the compound should have
            mol = self.to_pybel()
            self.rigid_inds = sorted(connect_rings(mol), key=lambda x: x[0])
        else:
            self.rigid_inds = None

        if self.name.endswith("typed"):
            # This is a hack to allow the old ff and typed files to work
            # We need to rename the atom types
            # TODO : test that gafffoyer and foyeroplsaa get same results and
            # remove this logic
            compound_pmd = pmd.load_file(input_str)
            for atom_pmd, atom_mb in zip(compound_pmd, self):
                atom_mb.name = "_{}".format(atom_pmd.type)


class Pack:
    def __init__(
        self,
        compound,
        n_compounds,
        density,
        ff=FORCE_FIELD["opv_gaff"],
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
        self.ff = ff
        self.remove_hydrogen_atoms = remove_hydrogen_atoms
        self.L = self._calculate_L()
        self.rigid_inds = []
        self.rigid_typeids = []

    def _remove_hydrogen(self):
        # TODO - not implemented with rigid
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

        L = self.L * box_expand_factor
        # Extra factor to make packing faster, will shrink it out
        box = mb.Box([L, L, L])
        system = mb.packing.fill_box(
            self.compound,
            n_compounds=self.n_compounds,
            box=box,
            overlap=0.2,
            fix_orientation=True,
        )

        # Calculate the rigid_inds in the packed system
        if any([comp.rigid_inds for comp in self.compound]):
            particle_count = 0
            rigid_count = 0
            for comp,n in zip(self.compound, self.n_compounds):
                if comp.rigid_inds is not None:
                    for _ in range(n):
                        for i,rigid in enumerate(comp.rigid_inds):
                            self.rigid_inds.append(rigid+particle_count)
                            self.rigid_typeids.append(i+rigid_count)
                        particle_count += comp.n_particles
                    rigid_count += len(comp.rigid_inds)
                else:
                    particle_count += n * comp.n_particles

        system.box = box
        pmd_system = system.to_parmed(residues=[self.residues])
        typed_system = self.ff.apply(
            pmd_system, assert_angle_params=False, assert_dihedral_params=False
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
