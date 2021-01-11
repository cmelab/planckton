import unyt as u

planckton_units = u.UnitSystem(
        'planckton', mass_unit='amu', length_unit='nm', time_unit='s'
        )
planckton_units['energy'] = 'kJ'

constants= {
        "avogadro": 6.022140857e23 / u.mol,
        "boltzmann": 1.38064852e-23 * u.J / u.K,
        }
