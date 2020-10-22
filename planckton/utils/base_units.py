import unyt as u

planckton_units = u.UnitSystem(
        'planckton', mass_unit='amu', length_unit='nm', time_unit='s'
        )
planckton_units['energy'] = 'kJ'

# base units are for sulfur
# TODO should we use autoscale option instead??
base_units = {
        "avogadro": 6.022140857e23 / u.mol,
        "boltzmann": 1.38064852e-23 * u.J / u.K,
        "mass": (32.06 * u.amu).in_base('planckton'),
        "energy": (1.046 * u.kJ / u.mol).in_base('planckton'),
        "length": (0.35635948725613575 * u.nm).in_base('planckton'),
        }
