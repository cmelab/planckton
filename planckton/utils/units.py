"""Utility functions for handling units."""

import unyt as u

planckton_units = u.UnitSystem(
    "planckton", mass_unit="amu", length_unit="nm", time_unit="s"
)
planckton_units["energy"] = "kJ"

constants = {
    "avogadro": 6.022140857e23 / u.mol,
    "boltzmann": 1.38064852e-23 * u.J / u.K,
}


def quantity_to_string(quantity):
    """Break a unyt.quantity into a string.

    Convert a unyt.quantity into a string containing its value and units in
    string format separated by "_". Useful for serialization. Division symbols,
    "/", which interfere with signac's workspace directory structure will be
    replaced with a dash, "-".

    IMPORTANT: This function expects one quantity, not an array.

    Parameters
    ----------
    quantity: unyt.unyt_quantity

    Returns
    -------
    str :
        unyt quantity in form "number_unit"
    """
    return f"{quantity.item()}_{str(quantity.units).replace('/','-')}"


def string_to_quantity(string):
    """Convert a string to unyt.quantity.

    Convert a string of  values and units separated by an underscore, "_", into
    a unyt quantity. See also `quantity_to_string`.

    Parameters
    ----------
    string : str
        Unyt quantity formatted as string, e.g., "number_units"

    Returns
    -------
    unyt.unyt_quantity
    """
    num, unit = string.split("_")
    return float(num) * u.Unit(unit.replace("-", "/"))


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
    T_SI = (
        T_reduced
        * ref_energy
        / (constants["boltzmann"] * constants["avogadro"])
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
        ref_mass * ref_distance ** 2 * constants["avogadro"] / ref_energy
    )
    real_time = dt * (time_squared ** 0.5)
    return real_time
