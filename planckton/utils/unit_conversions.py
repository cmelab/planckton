from planckton.utils.base_units import constants


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
