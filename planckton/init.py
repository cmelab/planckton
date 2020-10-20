import numpy as np

import mbuild as mb
import parmed as pmd
from planckton.utils import base_units


class Compound(mb.Compound):
    """ Wrapper class for mb.Compound"""

    def __init__(self, path_to_mol2):
        super(Compound, self).__init__()
        mb.load(path_to_mol2, compound=self)
        # Calculate mass of compound
        self.mass = np.sum([atom.mass for atom in self.to_parmed().atoms])
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
        out_file="init.hoomdxml",
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

        self.density = density
        self.ff_file = ff_file
        self.out_file = out_file
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
        Optional:
            box_expand_factor - float, Default = 5
            Expand the box before packing for faster
            packing.
        """
        units = base_units.base_units()

        if self.remove_hydrogen_atoms:
            self._remove_hydrogen()

        L = (
            self.L * box_expand_factor
        )  # Extra factor to make packing faster, will shrink it out
        box = mb.packing.fill_box(
            self.compound,
            n_compounds=self.n_compounds,
            # box=[-L/2, -L/2, -L/2, L/2, L/2, L/2],
            box=[L, L, L],
            overlap=0.2,
            edge=0.5,
            fix_orientation=True,
        )
        box.save(
            self.out_file,
            overwrite=True,
            forcefield_files=self.ff_file,
            ref_mass=units["mass"],  # amu
            ref_energy=units["energy"],  # kJ/mol
            ref_distance=units["distance"],  # nm
        )

    def _calculate_L(self):
        total_mass = np.sum(
            [n * c.mass for c, n in zip(self.compound, self.n_compounds)]
        )
        # Conversion from (amu/(g/cm^3))**(1/3) to ang
        L = (total_mass / self.density) ** (1 / 3) * 1.1841763
        L /= 10  # convert ang to nm
        return L


def test_typing(compound_file, ff_file):
    compound_mb = mb.load(compound_file)
    compound_pmd = pmd.load_file(compound_file)
    types_needed = set()
    for atom_pmd, atom_mb in zip(compound_pmd, compound_mb):
        atom_mb.name = "_{}".format(atom_pmd.type)
        types_needed.add(atom_pmd.type)
    box = mb.packing.fill_box(
        compound_mb,
        n_compounds=50,
        box=[6.53, 6.53, 6.53],
        overlap=0.2,
        fix_orientation=False,
    )
    box.save(
        "test_typing.hoomdxml",
        overwrite=True,
        forcefield_files="compounds/gaff.4fxml",
        ref_mass=32.06,
        ref_energy=1.046,
        ref_distance=0.35635948725613575,
    )
