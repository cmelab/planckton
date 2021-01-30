import unyt as u

from planckton.utils import units

ref_energy = 1.046 * u.kJ / u.mol
ref_mass = 32.06 * u.amu
ref_distance = 0.35635948725613575 * u.nm


def test_conversions():
    quantity = 1 * u.gram
    new_quantity = units.tuple_to_quantity(units.quantity_to_tuple(quantity))
    assert quantity == new_quantity

    tup = (1, "g")
    new_tup = units.quantity_to_tuple(units.tuple_to_quantity(tup))
    assert tup == new_tup

def test_temperature_reduction():
    # Compare conversion for SI (kelvin) to reduced kT
    T_reduced_true = 2.17
    T_SI = 273 * u.kelvin
    T_reduced = units.reduced_from_kelvin(T_SI, ref_energy)
    print(T_reduced)
    reduced_error = abs(T_reduced_true - T_reduced)
    assert reduced_error <= 0.001, "The error in the reduced kT is too high!"


def test_temperature_to_SI_conversion():
    # Compare conversion for kT (reduced) to SI (kelvin)
    T_SI_true = 126 * u.kelvin
    T_reduced = 1
    T_SI = units.kelvin_from_reduced(T_reduced, ref_energy)
    print(T_SI)
    SI_error = abs(T_SI_true - T_SI)
    assert SI_error <= 0.22, "The error in T (SI) is too high!"


def test_dt_to_SI_conversion():
    # Check real time
    timestep_true = 1.973 * u.fs
    dt_reduced = 0.001
    dt = units.convert_to_real_time(
            dt_reduced, ref_mass, ref_distance, ref_energy
            )
    print(dt.to('fs'))
    time_error = abs(timestep_true - dt)
    assert time_error <= 0.001, "The error in the timestep is too high!"


if __name__ == "__main__":
    test_temperature_reduction()
    test_temperature_to_SI_conversion()
    test_dt_to_SI_conversion()
