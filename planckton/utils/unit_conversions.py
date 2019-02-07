from planckton.utils import base_units


def reduce_from_kelvin(T_SI, precision=2):
    units = base_units.base_units()
    T = (units["avogadro"] * units["boltzmann"] * T_SI) / (
        units["kj_to_j"] * units["energy"]
    )
    T = round(T, precision)
    return T


def kelvin_from_reduced(T_reduced, precision=0):
    units = base_units.base_units()
    T_SI = (T_reduced * units["energy"] * units["kj_to_j"]) / (
        units["avogadro"] * units["boltzmann"]
    )
    T_SI = round(T_SI, precision)
    return T_SI


def convert_to_real_time(dt, precision=3):
    units = base_units.base_units()
    time_tau = units["mass"] * units["amu_to_kg"]
    time_tau *= (units["distance"] * units["nm_to_m"]) ** 2
    time_tau /= units["energy"] * units["kj_to_j"] / units["avogadro"]
    real_time = dt * (time_tau ** 0.5) * 1e15
    real_time = round(real_time, precision)
    return real_time
