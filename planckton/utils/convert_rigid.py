import json
import os
from collections import Counter

import hoomd
import hoomd.data
import hoomd.md
from hoomd.deprecated.init import read_xml
import numpy as np
from scipy.spatial import KDTree


def relative_pbc(positions, origin, box):
    box = np.array(box)
    p = np.copy(positions)
    p -= origin
    return shift_pbc(p,box)[0]

# TODO: I think hoomd can do this:
# http://hoomd-blue.readthedocs.io/en/stable/module-hoomd-data.html?highlight=wrap#hoomd.data.boxdim.min_image
def pbc_min_image(p1, p2, axes):
    dr = p1 - p2
    for i, p in enumerate(dr):
        if abs(dr[i]) > axes[i] * 0.5:
            p2[i] = (p2[i] + np.sign(dr[i]) * axes[i])
            # Use dr to decide if we need to add or subtract axis
    return p2


def com(points, masses):
    weighted = masses[:, None] * points
    M = np.sum(masses)
    return np.sum(weighted, axis=0) / M, M


def pbc_traslate(points, axes):
    """Translates a group of points to the minimum image of first point in list
    points: list of poits
    axes: array, Lx, Ly, Lz
    We assume poits are -L/2 to L/2, and a tetragonal unit cell
    By default, all points will be translated into the minimum image
    of the first point, but any point can be used
    """
    ref_point = points[0]
    # Add our min vector to our point
    min_image_cords = [
            pbc_min_image(ref_point, point, axes) for point in points[1:]
            ]
    # Skip over the first point since its our ref point
    min_image_cords.insert(0, ref_point)
    # Don't forget to add the ref point into the list of points
    return min_image_cords


def moit(points, masses):
    # TODO: Currently assumes center of mass is at origin, which is not
    # currently checked, so moments will be wrong.
    """Moment of Inertia Tensor
    Assumes center of mass is at origin
    Assumes 3xN array for points, 1xN for masses"""
    m = masses  # Makes things look nicer
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    I_xx = np.sum((y ** 2 + z ** 2) * m)
    I_yy = np.sum((x ** 2 + z ** 2) * m)
    I_zz = np.sum((x ** 2 + y ** 2) * m)
    return np.array((I_xx, I_yy, I_zz))


class RigidBody:
    pass


class NonRigidElements:
    pass


def set_topology(idx_map, old_snap, system):
    # Bonds
    new_bonds = [[idx_map[a], idx_map[b]] for (a, b) in old_snap.bonds.group]
    bond_types = {
            b_id: b_type for (b_id, b_type) in enumerate(old_snap.bonds.types)
            }
    for bond, bond_id in zip(new_bonds, old_snap.bonds.typeid):
        system.bonds.add(bond_types[bond_id], bond[0], bond[1])

    # Angles
    new_angles = [
        [idx_map[a], idx_map[b], idx_map[c]]
        for (a, b, c) in old_snap.angles.group
    ]
    angles_types = {
            a_id: a_type for (a_id, a_type) in enumerate(old_snap.angles.types)
            }
    for angle, angle_id in zip(new_angles, old_snap.angles.typeid):
        system.angles.add(angles_types[angle_id], angle[0], angle[1], angle[2])

    # Dihedrals
    new_dihedrals = [
        [idx_map[a], idx_map[b], idx_map[c], idx_map[d]]
        for (a, b, c, d) in old_snap.dihedrals.group
    ]
    dihedrals_types = {
        d_id: d_type for (d_id, d_type) in enumerate(old_snap.dihedrals.types)
    }
    for dihedral, dihedral_id in zip(new_dihedrals, old_snap.dihedrals.typeid):
        system.dihedrals.add(
            dihedrals_types[dihedral_id],
            dihedral[0],
            dihedral[1],
            dihedral[2],
            dihedral[3],
        )

    # Impropers
    new_impropers = [
        [idx_map[a], idx_map[b], idx_map[c], idx_map[d]]
        for (a, b, c, d) in old_snap.impropers.group
    ]
    impropers_types = {
        i_id: i_type for (i_id, i_type) in enumerate(old_snap.impropers.types)
    }
    for improper, improper_id in zip(new_impropers, old_snap.impropers.typeid):
        system.impropers.add(
            impropers_types[improper_id],
            improper[0],
            improper[1],
            improper[2],
            improper[3],
        )

    return system


def create_map(old_xyz, new_xyz, N_r_bodies):

    my_tree = KDTree(np.vstack(old_xyz))
    distance, old_index = my_tree.query(np.vstack(new_xyz[N_r_bodies:]))
    # Need to slice out rigid body centers
    assert [k for k, v in Counter(old_index).items() if v > 1] == []
    # Make sure there are no duplicates
    old_to_new = {
            key: value + N_r_bodies for (value, key) in enumerate(old_index)
            }
    return old_to_new


def init_wrapper(
    xmlfile,
    restart_rigid=False,
    rigid_flex_xyz_file="rigid_center_flex.xml",
    rigid_json_file="rigid_info.json",
):

    out_path, f_name = os.path.split(os.path.abspath(xmlfile))
    if restart_rigid:
        system = continue_rigid(
                xmlfile, out_path, rigid_flex_xyz_file, rigid_json_file
                )
        return system

    system = read_xml(filename=xmlfile, wrap_coordinates=True)
    system = new_rigid(system, out_path)
    return system


def continue_rigid(top_file, path, rigid_flex_xyz_file, rigid_json_file):
    """
    top_file: str
        hoomdxml file with all bond data
    rigid_flex_zyz_file: str
        hoomdxml with rigid centers and flex bodies, no topo data
    rigid_json_file: str
        json file with rigid body info
    """

    hoomd.context.initialize()
    system_old = read_xml(top_file, wrap_coordinates=True)
    snap_old = system_old.take_snapshot(all=True)
    xyz_old = snap_old.particles.position[:]
    r_bodies = set(snap_old.particles.body)
    # 4294967295 = -1 32-bit int
    r_bodies.discard(4294967295)
    N_r_bodies = len(r_bodies)

    hoomd.context.initialize()
    system = read_xml(path + "/" + rigid_flex_xyz_file, wrap_coordinates=True)
    # Add missing particle types
    missing_types = list(
            set(snap_old.particles.types) - set(system.particles.types)
            )
    [system.particles.types.add(missing_type) for missing_type in missing_types]
    # Add other missing data in this dumb way since you can't add new types with
    # system API
    snap = system.take_snapshot(all=True)
    snap.bonds.types = snap_old.bonds.types
    snap.angles.types = snap_old.angles.types
    snap.dihedrals.types = snap_old.dihedrals.types
    snap.impropers.types = snap_old.impropers.types
    system.restore_snapshot(snap)
    # Read in rigid body data
    with open(path + "/" + rigid_json_file, "r") as json_data:
        rigid_body_data = json.load(json_data)

    # Create rigid bodies
    rigid = hoomd.md.constrain.rigid()

    for rbody in rigid_body_data:
        rigid.set_param(
            str(rbody["r_type"]),
            positions=rbody["r_positions"],
            types=rbody["r_types"]
        )

    rigid.create_bodies()

    # Need new snap with correct XYZ since we just made new rigid bodies
    snap = system.take_snapshot(all=True)
    xyz = snap.particles.position[:]

    # fix topo data
    idx_map = create_map(xyz_old, xyz, N_r_bodies)
    system = set_topology(idx_map, snap_old, system)
    return system


def new_rigid(system, path):
    old_sys_xyz = []
    new_sys_xyz = []
    axis = [system.box.Lx, system.box.Ly, system.box.Lz]
    snapshot = system.take_snapshot(all=True)
    snapshot_old = system.take_snapshot(all=True)
    nonrigid = hoomd.group.nonrigid()
    # r_bodies is a list of rigid bodies, 0, 1, 2,...
    r_bodies = set(snapshot.particles.body)
    # 4294967295 = -1 32-bit int
    r_bodies.discard(4294967295)
    r_body_com = np.zeros((len(r_bodies), 3))
    r_body_m = np.zeros((len(r_bodies)))
    r_body_type = [""] * len(r_bodies)
    # Get xyz of each rigid body COM + make type list
    r_body_list = [RigidBody() for i in range(len(r_bodies))]
    f_body_list = [NonRigidElements() for i in range(len(nonrigid))]

    for atom in system.particles:
        old_sys_xyz.append(atom.position)

    for idx, atom in enumerate(nonrigid):
        f_body_list[idx].type = atom.type
        f_body_list[idx].mass = atom.mass
        f_body_list[idx].position = atom.position

    for idx, body_num in enumerate(r_bodies):
        r_body_name = "_R" + str(body_num)
        r_body_type[body_num] = r_body_name
        system.particles.types.add(r_body_name)
        # Get index for each rigid body member
        rigid_body_index = np.where(snapshot.particles.body == body_num)[0]
        r_xyz = []
        r_mass = []
        r_con_types = []
        r_con_pos = []
        # get xyz and mass of each rigid body member
        for index in rigid_body_index:
            r_xyz.append(snapshot.particles.position[index])
            r_mass.append(snapshot.particles.mass[index])
            r_con_types.append(
                snapshot.particles.types[snapshot.particles.typeid[index]]
            )

        r_xyz = np.stack(r_xyz)
        r_mass = np.array(r_mass)
        min_image_cords = pbc_traslate(r_xyz, axis)
        r_body_com[body_num], r_body_m[body_num] = com(min_image_cords, r_mass)

        for index in rigid_body_index:
            position = system.particles[int(index)].position
            position = np.array([position])
            # Add another array layer to make relative_pbc happy
            orign = r_body_com[body_num]
            r_con_pos.append(relative_pbc(position, orign, axis)[0])
            # remove extra array layer to make everything else happy

        r_body_list[idx].positions = r_con_pos
        r_body_list[idx].types = r_con_types
        r_body_list[idx].type = r_body_name
        r_body_list[idx].com = r_body_com[body_num]
        r_body_list[idx].moit = moit(np.stack(r_con_pos), r_mass)
        r_body_list[idx].mass = r_body_m[body_num]

    # First lets check to see how many unique type lists their are

    unique_types = list(set(map(tuple, [my_RB.types for my_RB in r_body_list])))
    rigid_body_map = {
        rb_types: "_R" + str(i) for i, rb_type in enumerate(unique_types)
    }
    for my_RB in r_body_list:
        my_RB.type = rigid_body_map[tuple(my_RB.types)]

    r_body_types = list(rigid_body_map.values())

    # Check case where r_body_types might be the same if they have the same
    # body type list, if assert passes build r_body_info dict
    # absolute(a - b) <= (atol + rtol * absolute(b))
    ATOL = 5e-5
    RTOL = 1e-4
    rigid_bodies_info = []
    for rbt in r_body_types:
        same_rbt = [my_RB for my_RB in r_body_list if my_RB.type == rbt]
        for i, my_RB in enumerate(same_rbt):
            assert np.allclose(
                np.array(same_rbt[0].positions).flatten(),
                np.array(my_RB.positions).flatten(),
                atol=ATOL,
                rtol=RTOL,
            )
        r_body_info = {}
        r_body_info["r_type"] = same_rbt[0].type
        r_body_info["r_types"] = same_rbt[0].types
        # Sorry for having a list of numpy arrays
        r_body_info["r_positions"] = np.array(same_rbt[0].positions).tolist()
        rigid_bodies_info.append(r_body_info)

    # Now we need to save the information for the next run
    with open(path + "/rigid_info.json", "w") as outfile:
        json.dump(rigid_bodies_info, outfile, indent=4)

    # Make snapshot of new system
    new_system = hoomd.data.make_snapshot(
        N=len(nonrigid) + len(r_bodies),
        box=snapshot.box,
        particle_types=snapshot.particles.types + r_body_types,
        bond_types=snapshot.bonds.types,
        angle_types=snapshot.angles.types,
        dihedral_types=snapshot.dihedrals.types,
        improper_types=snapshot.impropers.types,
        pair_types=snapshot.pairs.types,
        dtype="float",
    )

    hoomd.context.initialize()
    system = hoomd.init.read_snapshot(new_system)
    rigid = hoomd.md.constrain.rigid()

    for i, my_RB in enumerate(r_body_list):
        system.particles[i].type = my_RB.type
        system.particles[i].position = my_RB.com
        system.particles[i].moment_inertia = my_RB.moit
        system.particles[i].mass = my_RB.mass
        rigid.set_param(
                my_RB.type, positions=my_RB.positions, types=my_RB.types
                )
    for i, my_f in enumerate(f_body_list, start=len(r_bodies)):
        system.particles[i].type = my_f.type
        system.particles[i].position = my_f.position
        system.particles[i].mass = my_f.mass

    rigid.create_bodies()

    for atom in system.particles:
        new_sys_xyz.append(atom.position)
    # Now we need to figure out bond, angle, dihedrial
    # Create map
    old_to_new = create_map(old_sys_xyz, new_sys_xyz, len(r_bodies))

    # Use map to set bonds, angles, dihedrals, and impropers
    system = set_topology(old_to_new, snapshot_old, system)
    return system
