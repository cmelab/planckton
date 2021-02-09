import unyt as u


planckton_units = u.UnitSystem(
        'planckton', mass_unit='amu', length_unit='nm', time_unit='s'
        )
planckton_units['energy'] = 'kJ'

constants= {
        "avogadro": 6.022140857e23 / u.mol,
        "boltzmann": 1.38064852e-23 * u.J / u.K,
        }


def quantity_to_tuple(quantity):
    """
    Break a unyt.quantity into value and units in string format.

    IMPORTANT: This function expects one quantity, not an array.

    Parameters
    ----------
    quantity: unyt.unyt_quantity

    Returns
    -------
    (number, string)
    """
    return (quantity.item(), str(quantity.units))

def tuple_to_quantity(tup):
    """
    Convert a tuple containing values and units in string format into a unyt
    quantity.

    Parameters
    ----------
    tup: tuple
        first value is a number and the second value is a string with the units

    Returns
    -------
    unyt.unyt_quantity
    """
    return tup[0] * u.Unit(tup[1])


def reduced_from_kelvin(T_SI, ref_energy):
    T = constants["boltzmann"] * constants["avogadro"] * T_SI / ref_energy
    return T


def kelvin_from_reduced(T_reduced, ref_energy):
    T_SI = T_reduced * ref_energy / (
            constants["boltzmann"] * constants["avogadro"]
            )
    return T_SI


def convert_to_real_time(dt, ref_mass, ref_distance, ref_energy):
    time_squared = (
            ref_mass * ref_distance**2 * constants["avogadro"] / ref_energy
            )
    real_time = dt * (time_squared ** 0.5)
    return real_time
