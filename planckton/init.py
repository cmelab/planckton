"""Tools for initializing an OPV simulation with PlanckTon."""

import os
import warnings

import foyer
import mbuild as mb
import numpy as np
import parmed as pmd
import unyt as u
from ele import element_from_symbol
from ele.exceptions import ElementError
from mbuild.coordinate_transform import z_axis_transform
from unyt.exceptions import UnitConversionError

from planckton.forcefields import FORCEFIELD
from planckton.utils.units import planckton_units


class Compound(mb.Compound):
    """Wrapper class for mb.Compound.

    Parameters
    ----------
    input_str : str
        Path to the file to load or SMILES string

    Attributes
    ----------
    mass : unyt.unyt_quantity
        The mass of the compound in amus
    name : str
        Compound name, used to apply the forcefield more quickly in foyer
    typed : bool
        Whether the compound is already typed
    """

    def __init__(self, input_str):
        super(Compound, self).__init__()
        if os.path.exists(input_str):
            mb.load(input_str, compound=self)
        else:
            mb.load(input_str, smiles=True, compound=self)

        # Calculate mass of compound
        self.set_elements()

        # This helps to_parmed use residues to apply ff more quickly
        self.name = os.path.basename(input_str).split(".")[0]

        if self.name.endswith("typed"):
            self.typed = True
            # This is a hack to allow the old ff and typed files to work
            # We need to rename the atom types
            # TODO : test that gafffoyer and foyeroplsaa get same results and
            # remove this logic
            compound_pmd = pmd.load_file(input_str)
            for atom_pmd, atom_mb in zip(compound_pmd, self):
                atom_mb.name = "_{}".format(atom_pmd.type)
        else:
            self.typed = False

    def set_elements(self):
        """Set the ele.element of each particle in the compound."""
        for p in self.particles():
            try:
                p.element = element_from_symbol(p.name)
            except ElementError:
                # This is a hack for our typed mol2 files
                p.element = element_from_symbol(p.name.strip("0123456789"))


class Pack:
    """Convenience class for filling box and atomtyping.

    Parameters
    ----------
    compound : Compound or list of Compounds
        Compound(s) to initialize in simulation
    n_compounds : int or list of ints
        Number(s) of compound(s) to initialize
    density : unyt.unyt_quantity
        Density of the system with units::
            import unyt as u
            1.0 * u.g / u.cm**3
    ff : foyer.Forcefield
        Foyer forcefield to use for typing compounds
        (default foyer.Forcefield("gaff-custom.xml"))
    remove_hydrogen_atoms : bool
        Whether to remove hydrogen atoms. (default False)
    foyer_kwargs = dict
        Keyword arguments to be passed to foyer.Forcefield.apply()
        (default {"assert_dihedral_params": False})

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
        ff=FORCEFIELD["gaff-custom"],
        remove_hydrogen_atoms=False,
        foyer_kwargs={"assert_dihedral_params": False},
    ):
        if not isinstance(compound, (list, set)):
            self.compound = [compound]
        else:
            self.compound = compound
        if n_compounds is not None and not isinstance(n_compounds, (list, set)):
            self.n_compounds = [n_compounds]
        else:
            self.n_compounds = n_compounds

        if not all([c.typed for c in self.compound]) and ff in builtin_ffs:
            raise NotImplementedError(
                "Untyped compounds not supported with custom forcefields. "
                "Try a forcefield with SMARTS grammar, like 'gaff', instead."
            )
        if any([c.typed for c in self.compound]) and ff not in builtin_ffs:
            warnings.warn(
                "You are using a typed compound. Make sure that the types "
                "match those in your forcefield. If you run into errors and "
                "are using 'gaff', try 'gaff-custom' instead."
            )

        if isinstance(density, u.unyt_quantity):
            try:
                # catch unit errors early
                density.to(
                    planckton_units["mass"] / planckton_units["length"] ** 3
                )
            except UnitConversionError as e:
                raise (e)
            self.density = density
        else:
            raise TypeError("density must be a unyt quantity")

        self.residues = [comp.name for comp in self.compound]
        self.ff = ff
        self.remove_hydrogen_atoms = remove_hydrogen_atoms
        self.L = self._calculate_L()
        self.foyer_kwargs = foyer_kwargs

    def pack(self, box_expand_factor=5):
        """Pack compounds into a larger box in preparation for shrinking.

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
        L = self.L.value * box_expand_factor
        box = mb.Box([L, L, L])
        system = mb.packing.fill_box(
            self.compound,
            n_compounds=self.n_compounds,
            box=box,
            overlap=0.2,
            fix_orientation=True,
        )
        system.box = box
        pmd_system = system.to_parmed(residues=[self.residues])
        typed_system = self.ff.apply(pmd_system, **self.foyer_kwargs)
        if self.remove_hydrogen_atoms:
            typed_system.strip([a.atomic_number == 1 for a in pmd_system.atoms])
        return typed_system

    def stack(self, separation=0.7):
        """Stack compounds in a single layer.

        Parameters
        ----------
        separation : float, default=0.7
            Distance in nm to separate the individually stacked compounds

        Returns
        -------
        typed_system : ParmEd structure
            ParmEd structure of the box with the stacked compounds
        """
        system = mb.Compound()
        L = self.L.value

        for idx, comp in enumerate(self.compound):
            z_axis_transform(comp)
            comp.translate(np.array([separation, 0, 0]) * idx)
            system.add(comp)

        system.box = mb.box.Box([L, L, L])

        system.translate_to(
            (system.box.Lx / 2, system.box.Ly / 2, system.box.Lz / 2)
        )

        pmd_system = system.to_parmed(residues=[self.residues])
        typed_system = self.ff.apply(pmd_system, **self.foyer_kwargs)
        if self.remove_hydrogen_atoms:
            typed_system.strip([a.atomic_number == 1 for a in pmd_system.atoms])
        return typed_system

    def _calculate_L(self):
        total_mass = (
            np.sum(
                [c.mass * n for c, n in zip(self.compound, self.n_compounds)]
            )
            * u.amu
        )

        L = (total_mass / self.density) ** (1 / 3)
        return L.in_base("planckton")


builtin_ffs = [FORCEFIELD["gaff-custom"], FORCEFIELD["oplsua-custom"]]
