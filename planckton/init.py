import numpy as np
import mbuild as mb
import parmed as pmd
import unyt as u
from unyt.exceptions import UnitConversionError

from planckton.utils.units import planckton_units


class Compound(mb.Compound):
    """
    Wrapper class for mb.Compound

    Parameters
    ----------
    path_to_mol2 : str
        Path to the mol2 file to load

    Attributes
    ----------
    mass : unyt.unyt_quantity
        The mass of the compound in amus
    """

    def __init__(self, path_to_mol2):
        super(Compound, self).__init__()
        mb.load(path_to_mol2, compound=self)
        # Calculate mass of compound
        # ParmEd uses amu
        # TODO use ele?
        self.mass = np.sum([a.mass for a in self.to_parmed().atoms]) * u.amu
        # We need to rename the atom types
        compound_pmd = pmd.load_file(path_to_mol2)
        for atom_pmd, atom_mb in zip(compound_pmd, self):
            atom_mb.name = "_{}".format(atom_pmd.type)


class Pack:
    """
    Convenience class for filling box and atomtyping

    Parameters
    ----------
    compound : Compound or list of Compounds
        Compound(s) to initialize in simulation
    n_compounds : int or list of ints
        Number(s) of compound(s) to initialize
    density : float or unyt.unyt_quantity
        Density of the system. If only a float is provided, the density is
        assumed to be in planckton units (amu/nm^3).
    ff_file : str
        Foyer forcefield xml file to use for typing compounds
        (default "compounds/gaff.4fxml")
    out_file : str
        Filename to write out the typed system (default "init.hoomdxml")
    remove_hydrogen_atoms : bool
        Whether to remove hydrogen atoms. (default False)

    Attributes
    ----------
    compound : list of Compound(s)
        Compound(s) in the system
    n_compounds : list of int(s)
        Number(s) of Compound(s) in the system
    density : unyt.unyt_quantity
        Density of the system with units
    ff_file : str
        Path to foyer forcefield xml
    out_file : str
        Path to where the hoomdxml where the initialized system will be written
    remove_hydrogen_atoms : bool
        Whether hydrogen atoms will be removed.
    L : unyt.unyt_quantity
        Length of the box with units
    """
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

        if isinstance(density, u.unyt_quantity):
            try:
                # catch unit errors early
                density.to(
                        planckton_units["mass"] / planckton_units["length"]**3
                        )
            except UnitConversionError as e:
                raise(e)
            self.density = density
        else:
            self.density = (
                    density *
                    planckton_units["mass"] / planckton_units["length"]**3
                    )

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

        if self.remove_hydrogen_atoms:
            self._remove_hydrogen()

        L = (self.L.value * box_expand_factor)
        box = mb.packing.fill_box(
            self.compound,
            n_compounds=self.n_compounds,
            box=[L, L, L],
            overlap=0.2,
            edge=0.5,
            fix_orientation=True,
        )
        box.save(
            self.out_file,
            overwrite=True,
            forcefield_files=self.ff_file,
            auto_scale=True,
            foyer_kwargs={"assert_dihedral_params": False}
        )

    def _calculate_L(self):
        total_mass = np.sum([
            n * c.mass.in_base('planckton')
            for c, n in zip(self.compound, self.n_compounds)
            ]) * u.amu

        L = (total_mass / self.density) ** (1 / 3)
        return L.in_base('planckton')

