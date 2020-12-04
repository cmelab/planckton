import hoomd
from hoomd.data import make_snapshot
from mbuild.formats.hoomd_simulation import create_hoomd_simulation
import numpy as np


def init_rigid(rigid_inds, rigid_typeids, typed_system, sim):
    # Determine how many rigid bodies the system should have
    n_bodies = len(rigid_inds)
    with sim:
        # Make an initial snapshot with only general rigid body centers
        # --the box length doesn't matter and specific types are added after
        init_snap = make_snapshot(
            N=n_bodies, particle_types=["_R"], box=hoomd.data.boxdim(L=10)
        )
        # set the rigid types and typeids
        init_snap.particles.types = [
                f"_R{i}" for i in range(max(rigid_typeids)+1)
                ]
        init_snap.particles.typeid[:] = rigid_typeids

        # Add the typed system to this snapshot
        hoomd_objects, ref_values = create_hoomd_simulation(
            typed_system, auto_scale=True, init_snap=init_snap
        )
        snap = hoomd_objects[0]

        for i, ring in enumerate(rigid_inds):
            # Indices of constituent particles
            inds = ring + n_bodies

            # Move the rigid body centers
            snap.particles.position[i] = np.mean(
                snap.particles.position[inds], axis=0
            )

            # Set body tags
            snap.particles.body[i] = i
            snap.particles.body[inds] = i * np.ones(len(ring))

            # Moment of inertia tensor
            snap.particles.moment_inertia[i] = moit(
                snap.particles.position[inds],
                snap.particles.mass[inds],
                center=snap.particles.position[i],
            )

            # Mass
            snap.particles.mass[i] = np.sum(snap.particles.mass[inds])

        # Reinitialize with modified snapshot
        sim.system_definition.initializeFromSnapshot(snap)

        # Add body exclusions to neighborlist
        nl = sim.neighbor_lists[0]
        ex_list = nl.exclusions
        ex_list.append("body")
        sim.neighbor_lists[0].reset_exclusions(exclusions=ex_list)

        # Set the body types for each rigid body type
        rigid = hoomd.md.constrain.rigid()
        for i in range(max(rigid_typeids)+1):
            # particles indices of the first instance center and constituent
            # body for this rigid type
            first_center = np.where(snap.particles.typeid == i)[0][0]
            first_body = rigid_inds[rigid_typeids.index(i)]

            center_pos = snap.particles.position[first_center]
            # body indices must be shifted by the added center particles
            body_pos = snap.particles.position[first_body + n_bodies]
            body_pos -= center_pos

            body_types = [
                    snap.particles.types[i]
                    for i in snap.particles.typeid[first_body + n_bodies]
                    ]
            rigid.set_param(
                    f"_R{i}",
                    types=body_types,
                    positions=[tuple(i) for i in body_pos]
                    )
        rigid.validate_bodies()

        # add zero interactions for rigid body centers with all particles
        lj = sim.forces[0]
        lj.pair_coeff.set(
                [f"_R{i}" for i in range(max(rigid_typeids)+1)],
                snap.particles.types,
                epsilon=0,
                sigma=0
                )

        centers = hoomd.group.rigid_center()
        nonrigid = hoomd.group.nonrigid()
        both = hoomd.group.union("both", centers, nonrigid)
    return both, snap, sim, ref_values


def connect_rings(mol):
    """
    Get the particle indices of the conjugated particles in the system.
    This is useful because conjugated systems are expected to remain planar
    --these indices can be used to specify the rigid bodies.

    Parameters
    ----------
    mol : openbabel.pybel.Molecule
        A molecule object containing conjugated rings

    Returns
    -------
    ring_arrs : list of numpy.arrays
        Each array contains the particle indices of each connected
        conjugated system.
    """
    # SSSR - find Smallest Set of Smallest Rings
    # http://openbabel.org/dev-api/classOpenBabel_1_1OBRing.shtml
    # pybel indices start at 1, so they are shifted to match every
    # other python library
    rings = [set(np.array(ring._path) - 1) for ring in mol.OBMol.GetSSSR()]

    # Iterate through rings until they are all connected
    connected = False
    while not connected:
        rings, connected = _check_rings(rings)

    # convert to numpy array so it can be used for indexing
    ring_arrs = []
    for ring in rings:
        ring_arrs.append(np.array([*ring]))
    return ring_arrs


def _check_rings(rings):
    # if not all rings are disjoint, then some must still share particles
    connected = all(
        [
            ringi.isdisjoint(ringj)
            for i, ringi in enumerate(rings[:-1])
            for ringj in rings[i + 1 :]
        ]
    )
    if not connected:
        new_rings = [
            set(sorted(ringi.union(ringj)))
            for i, ringi in enumerate(rings[:-1])
            for ringj in rings[i + 1 :]
            if not ringi.isdisjoint(ringj)
        ]

        # this ends up adding each connected ring twice, so the next
        # section fixes that
        conjugated = []
        for i in new_rings:
            if i not in conjugated:
                conjugated.append(i)

        # Add any disjoint rings that are already fully connected
        for i in rings:
            disjoint = [i.isdisjoint(j) for j in conjugated]
            if all(disjoint):
                if sorted(i) not in conjugated:
                    conjugated.append(set(sorted(i)))
    else:
        conjugated = rings
    return conjugated, connected


def moit(points, masses, center=np.zeros(3)):
    """
    Calculates moment of inertia tensor (moit) for rigid bodies. Assumes
    rigid body center is at origin unless center is provided.
    Only calculates diagonal elements.

    Parameters
    ----------
    points : numpy.ndarray (N,3)
        x, y, and z coordinates of the rigid body constituent particles
    masses : numpy.ndarray (N,)
        masses of the constituent particles
    center : numpy.ndarray (3,)
        x, y, and z coordinates of the rigid body center
        (default np.array([0,0,0]))

    Returns
    -------
    numpy.ndarray (3,)
        moment of inertia tensor for the rigid body center
    """
    points -= center
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    I_xx = np.sum((y ** 2 + z ** 2) * masses)
    I_yy = np.sum((x ** 2 + z ** 2) * masses)
    I_zz = np.sum((x ** 2 + y ** 2) * masses)
    return np.array((I_xx, I_yy, I_zz))
