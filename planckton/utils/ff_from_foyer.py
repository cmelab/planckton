import xml.etree.cElementTree as ET
from collections import OrderedDict
from itertools import tee

import numpy as np

import hoomd
import hoomd.md


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def parse_more_than_one_period(periods):
    flat_list = []
    # dihedral type
    flat_list.append(periods[0][0])
    # build flat list with all periods
    for coeffs in periods:
        # slice off repeated type
        for parm in coeffs[1:]:
            flat_list.append(parm)
    return flat_list


def periodic_torsion_force(
    theta,
    psi_k1,
    period1,
    phase1,
    psi_k2=0,
    period2=0,
    phase2=0,
    psi_k3=0,
    period3=0,
    phase3=0,
    psi_k4=0,
    period4=0,
    phase4=0,
):
    V = (
        psi_k1 * (1 + np.cos((period1 * theta) + phase1))
        + psi_k2 * (1 + np.cos((period2 * theta) + phase2))
        + psi_k3 * (1 + np.cos((period3 * theta) + phase3))
        + psi_k4 * (1 + np.cos((period4 * theta) + phase4))
    )
    F = (
        psi_k1 * period1 * np.sin((period1 * theta) + phase1)
        + psi_k2 * period2 * np.sin((period2 * theta) + phase2)
        + psi_k3 * period3 * np.sin((period3 * theta) + phase3)
        + psi_k4 * period4 * np.sin((period4 * theta) + phase4)
    )
    return (V, F)


def set_coeffs(file_name, system, nl, e_factor=1.0, r_cut=2.5):
    """
    Read in the molecular dynamics coefficients exported by Foyer

    Parameters
    ----------
    file_name : str
        Path to foyer forcefield xml
    system : Hoomd system object
        The system to set coefficients against
    nl : Hoomd neighborlist object
        The neighbourlist
    """

    coeffs_dict = get_coeffs(file_name)
    lj = hoomd.md.pair.lj(r_cut=r_cut, nlist=nl)
    lj.set_params(mode="xplor")

    coeffs_dict["pair_coeffs"] += [
        [_, 0.0, 0.0] for _ in system.particles.types if _.startswith("_R")
    ]

    for type1 in coeffs_dict["pair_coeffs"]:
        for type2 in coeffs_dict["pair_coeffs"]:
            lj.pair_coeff.set(
                type1[0],
                type2[0],
                epsilon=e_factor * np.sqrt(type1[1] * type2[1]),
                sigma=np.sqrt(type1[2] * type2[2]),
            )

    if coeffs_dict["bond_coeffs"]:
        harmonic_bond = hoomd.md.bond.harmonic()
        for bond in coeffs_dict["bond_coeffs"]:
            harmonic_bond.bond_coeff.set(bond[0], k=bond[1], r0=bond[2])

    if coeffs_dict["angle_coeffs"]:
        harmonic_angle = hoomd.md.angle.harmonic()
        for angle in coeffs_dict["angle_coeffs"]:
            harmonic_angle.angle_coeff.set(angle[0], k=angle[1], t0=angle[2])

    if coeffs_dict["dihedral_coeffs"]:
        if len(coeffs_dict["dihedral_coeffs"][0]) == 5:  # OPLS Check
            harmonic_dihedral = hoomd.md.dihedral.opls()
            for dihedral in coeffs_dict["dihedral_coeffs"]:
                harmonic_dihedral.dihedral_coeff.set(
                    dihedral[0],
                    k1=dihedral[1],
                    k2=dihedral[2],
                    k3=dihedral[3],
                    k4=dihedral[4],
                )
        else:  # AMBER style
            dtable = hoomd.md.dihedral.table(width=1000)
            group = []
            for cur, _next in pairwise(coeffs_dict["dihedral_coeffs"]):
                if cur[0] == _next[0]:
                    group.append(cur)
                else:
                    group.append(cur)
                    if len(group) > 1:
                        parsed_di = parse_more_than_one_period(group)
                        my_coeffs = OrderedDict(
                            [
                                ("psi_k1", 0),
                                ("period1", 0),
                                ("phase1", 0),
                                ("psi_k2", 0),
                                ("period2", 0),
                                ("phase2", 0),
                                ("psi_k3", 0),
                                ("period3", 0),
                                ("phase3", 0),
                                ("psi_k4", 0),
                                ("period4", 0),
                                ("phase4", 0),
                            ]
                        )
                        for _ in range(len(parsed_di[1:]) // 3):
                            for key, val in zip(
                                my_coeffs, parsed_di[1:][_ * 3 : _ * 3 + 3]
                            ):
                                key = key[:-1] + str(int(key[-1]) + _)
                                my_coeffs[key] = val
                        dtable.dihedral_coeff.set(
                            parsed_di[0],
                            func=periodic_torsion_force,
                            coeff=my_coeffs
                        )
                    else:
                        # call set di here, only 1
                        my_coeffs = {
                            "psi_k1": cur[1],
                            "period1": cur[2],
                            "phase1": cur[3],
                        }
                        dtable.dihedral_coeff.set(
                            cur[0], func=periodic_torsion_force, coeff=my_coeffs
                        )
                    group = []
            else:
                group.append(cur)
                if len(group) > 1:
                    parse_more_than_one_period(group)
                else:
                    my_coefs = {
                        "psi_k1": _next[1],
                        "period1": _next[2],
                        "phase1": _next[3],
                    }
                    dtable.dihedral_coeff.set(
                        _next[0], func=periodic_torsion_force, coeff=my_coeffs
                    )

    for atomID, atom in enumerate(system.particles):
        if not str(atom.type).startswith("_R"):
            atom.mass = coeffs_dict["mass"][str(atom.type)]
    # TODO: Support for improppers
    # TODO: Support for charges
    # pppmnl = hoomd.md.nlist.cell()
    # pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    # pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)

    return system


def get_coeffs(file_name):
    coeff_dictionary = {
        "pair_coeffs": [],
        "bond_coeffs": [],
        "angle_coeffs": [],
        "dihedral_coeffs": [],
    }
    with open(file_name, "r") as xml_file:
        xml_data = ET.parse(xml_file)
    root = xml_data.getroot()
    for config in root:
        for child in config:
            # First get the masses which are different
            if child.tag == "mass":
                masses = [
                        float(_) for _ in child.text.split("\n") if len(_) > 0
                        ]
            # Now the other coefficients
            if child.tag == "type":
                types = [str(_) for _ in child.text.split("\n") if len(_) > 0]
            elif child.tag in coeff_dictionary.keys():
                if child.text is None:
                    continue
                for line in child.text.split("\n"):
                    if len(line) == 0:
                        continue
                    coeff = line.split()
                    coeff_dictionary[child.tag].append(
                        [coeff[0]] + list(map(float, coeff[1:]))
                    )
    all_mass_type_pairs = [[pair[0], pair[1]] for pair in zip(types, masses)]
    coeff_dictionary["mass"] = {
        unique[0]: unique[1]
        for unique in set(tuple(unique) for unique in all_mass_type_pairs)
    }
    return coeff_dictionary
