"""Utility functions for handling units."""

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
    """Break a unyt.quantity into a tuple.

    Convert a unyt.quantity into a tuple containing its value and units in
    string format. Useful for serialization.

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
    """Convert tuple to unyt.quantity.

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
    """Convert temperature in Kelvin to reduced temperature.

    Parameters
    ----------
    T_SI: unyt.quantity
        Temperature in Kelvin.
    ref_energy: unyt.quantity
        Reference energy from the simulation

    Returns
    -------
    unyt.unyt_quantity
        Unitless temperature
    """
    T = constants["boltzmann"] * constants["avogadro"] * T_SI / ref_energy
    return T


def kelvin_from_reduced(T_reduced, ref_energy):
    """Convert temperature in reduced units to Kelvin.

    Parameters
    ----------
    T_reduced: unyt.quantity or float
        Temperature in simulation units. (unitless)
    ref_energy: unyt.quantity
        Reference energy from the simulation

    Returns
    -------
    unyt.unyt_quantity
        Temperature in Kelvin
    """
    T_SI = T_reduced * ref_energy / (
            constants["boltzmann"] * constants["avogadro"]
            )
    return T_SI


def convert_to_real_time(dt, ref_mass, ref_distance, ref_energy):
    """Convert the timestep in reduced units to seconds.

    Parameters
    ----------
    dt: float
        Timestep
    ref_mass: unyt.quantity
        Reference mass from the simulation
    ref_distance: unyt.quantity
        Reference distance from the simulation
    ref_energy: unyt.quantity
        Reference energy from the simulation

    Returns
    -------
    unyt.unyt_quantity
        Timestep in seconds
    """
    time_squared = (
            ref_mass * ref_distance**2 * constants["avogadro"] / ref_energy
            )
    real_time = dt * (time_squared ** 0.5)
    return real_time
