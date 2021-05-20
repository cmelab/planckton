"""Utility functions for scaling energies for implicit solvents."""


def set_coeffs(pair_dict, e_factor):
    """Scale energy coefficients.

    Given an e_factor and a dictionary of pair parameters (e.g. from
    `lj.get_metadata()['pair_coeff'].get_metadata()` where `lj` is
    `hoomd.md.pair.lj`). This function will return parameters to pass into
    the set function for a pair potential (e.g. `lj.set(a, b, **coeffs)`).

    Parameters
    ----------
    pair_dict : dict
        output from `lj.get_metadata()['pair_coeff'].get_metadata()`
    e_factor : float
        Scaling factor multiplied to the epsilon value in the dictionary

    Returns
    -------
    a,b : str
        Type names for particles i and j
    new_dict : dict
        Copy of pair_dict with typei/typej keys removed and epsilon scaled
        by e_factor
    """
    new_dict = dict(pair_dict)
    a = new_dict.pop("typei")
    b = new_dict.pop("typej")
    new_dict["epsilon"] *= e_factor
    return a, b, new_dict
