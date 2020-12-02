import numpy as np


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
    rings = [set(np.array(ring._path)-1) for ring in mol.OBMol.GetSSSR()]

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
            for i,ringi in enumerate(rings[:-1])
            for ringj in rings[i+1:]
        ]
    )
    if not connected:
        conjugated = [
            set(sorted(ringi.union(ringj)))
            for i,ringi in enumerate(rings[:-1])
            for ringj in rings[i+1:]
            if not ringi.isdisjoint(ringj)
        ]

        # this ends up adding each connected ring twice, so the next
        # section fixes that
        res = []
        [res.append(i) for i in conjugated if i not in res]
    else:
        res = rings
    return res, connected

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
