from planckton.utils import unit_conversions


def check_temperature_reduction():
    # Compare conversion for SI (kelvin) to reduced kT
    T_reduced_true = 2.17
    T_SI = 273  # K
    T_reduced = unit_conversions.reduce_from_kelvin(T_SI)
    reduced_error = abs(T_reduced_true - T_reduced)
    assert reduced_error <= 0.001, "The error in the reduced kT is too high!"


def check_temperature_to_SI_conversion():
    # Compare conversion for kT (reduced) to SI (kelvin)
    T_SI_true = 126  # K
    T_reduced = 1
    T_SI = unit_conversions.kelvin_from_reduced(T_reduced)
    SI_error = abs(T_SI_true - T_SI)
    assert SI_error <= 0.001, "The error in T (SI) is too high!"


def check_dt_to_SI_conversion():
    # Check real time
    timestep_true = 1.973  # fs
    dt_reduced = 0.001
    dt = unit_conversions.convert_to_real_time(dt_reduced)
    time_error = abs(timestep_true - dt)
    assert time_error <= 0.001, "The error in the timestep is too high!"


if __name__ == "__main__":
    check_temperature_reduction()
    check_temperature_to_SI_conversion()
    check_dt_to_SI_conversion()
