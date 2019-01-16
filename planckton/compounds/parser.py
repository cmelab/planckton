from copy import deepcopy
from collections import defaultdict

PI = 3.141_592_653_589


class AmberFFParams:
    types = []
    bonds = []
    angles = []
    dihedrals = []
    impropers = []
    lj = []


def kcal2kJ(val):
    return float(val) * 4.184


def ang2nm(val):
    return float(val) / 10


def amber_bond2openmm_bond(val):
    # don't forget the factor of two between functional forms
    # 100 is for 1/ang to 1/nm
    return kcal2kJ(val) * (100 * 2)


def amber_angle2openmm_angle(val):
    # don't forget the factor of two between functional forms
    return kcal2kJ(val) * 2


def deg2rad(val):
    return float(val) * (PI / 180)


def parse_mass(line):
    line = line.strip().split(" ")
    line = list(filter(None, line))
    amber_type, mass, atomic_polarizability = line
    AtomType(amber_type, mass, atomic_polarizability)


def parse_bond(line):
    line = line.strip().split(" ")
    line = list(filter(None, line))
    if "-" not in line[0]:
        line[:2] = ["".join(line[:2])]
    classes = line[0].split("-")
    k, r0 = line[1:3]
    # Unit conversions
    k, r0 = amber_bond2openmm_bond(k), ang2nm(r0)
    HarmonicBondForce(*classes, k, r0)


def parse_angle(line):
    line = line.strip().split(" ")
    line = list(filter(None, line))
    while len(line[0].split("-")) != 3:
        line[0] += line.pop(1)
    classes = line[0].split("-")
    k, theta0 = line[1:3]
    # Unit conversions
    k, theta0 = amber_angle2openmm_angle(k), deg2rad(theta0)
    HarmonicAngleForce(*classes, k, theta0)


def parse_dihedral(line):
    line = line.strip().split(" ")
    line = list(filter(None, line))
    while len(line[0].split("-")) != 4:
        line[0] += line.pop(1)
    classes = line[0].split("-")
    devider, barrier, phase, periodicity = map(float, line[1:5])
    # http://ambermd.org/doc12/Amber18.pdf#page=247&zoom=auto,-398,716
    # Section 14.1.6 in amber18 manual
    k_eff = barrier / devider
    # Unit conversions
    k_eff, phase = kcal2kJ(k_eff), deg2rad(phase)
    PeriodicTorsionForce(*classes, k_eff, phase, periodicity)


def parse_improper(line):
    line = line.strip().split(" ")
    line = list(filter(None, line))
    while len(line[0].split("-")) != 4:
        line[0] += line.pop(1)
    classes = line[0].split("-")
    barrier, phase, periodicity = map(float, line[1:4])
    k_eff = barrier
    # Unit conversions
    k_eff, phase = kcal2kJ(k_eff), deg2rad(phase)
    PeriodicTorsionForceImproper(*classes, k_eff, phase, periodicity)


def parse_lj(line):
    line = line.strip().split(" ")
    line = list(filter(None, line))
    amber_type, r_min, epsilon = line
    # Unit conversions
    sigma, epsilon = ang2nm(r_min), kcal2kJ(epsilon)
    sigma = sigma * 2 ** (-1 / 6) * 2
    NonbondedForce(amber_type, sigma, epsilon)


sections = {
    "mass": {"key": "MASS\n", "parser": parse_mass},
    "bond": {"key": "BOND\n", "parser": parse_bond},
    "angle": {"key": "ANGLE\n", "parser": parse_angle},
    "dihedral": {"key": "DIHE\n", "parser": parse_dihedral},
    "improper": {"key": "IMPROPER\n", "parser": parse_improper},
    "lj": {"key": "NONBON\n", "parser": parse_lj},
}


class OpenMMXMLField:
    subclasses = []
    open_xml_prams = defaultdict(list)

    def gen_xml(self):
        pass

    def gen_xml_header(self):
        return f" <{self.xml_field}>\n"

    def gen_xml_closer(self):
        return f" </{self.xml_field}>\n"

    def __init__(self, *args, **kwargs):
        # Creation of an object that is a OpenMMXMLField or subclass of OpenMMXMLField
        if self.__class__ in (OpenMMXMLField,):
            # Ignore containers that are not real Pets
            pass
        else:
            if self not in self.open_xml_prams[self.__class__]:
                self.open_xml_prams[self.__class__].append(self)

    def __init_subclass__(cls, **kwargs):
        # Creation of a new class that inherits OpenMMXMLField
        super().__init_subclass__(**kwargs)
        cls.subclasses.append(cls)

    # TODO Strings are unquoted with this, need to fix
    def __repr__(self):
        attrs = ", ".join(f"{k} = {v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__name__}({attrs})"

    def __hash__(self):
        return hash(tuple(self.__dict__.values()) + (self.__class__.__name__,))

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False


class AtomType(OpenMMXMLField):
    xml_field = "AtomTypes"

    def __init__(self, name, mass, atomic_polarizability):
        self.name = name
        self.mass = mass
        self._class = name
        # self.element = name[0].upper()
        self.element = f"_{name}"
        self.atomic_polarizability = atomic_polarizability
        super().__init__()

    def gen_xml(self):
        return f'  <Type name="{self.name}" class="{self._class}" element="{self.element}" mass="{self.mass}" def="{self.element}" desc="" doi=""/>\n'


class HarmonicBondForce(OpenMMXMLField):
    xml_field = "HarmonicBondForce"

    def __init__(self, class_1, class_2, k, r0):
        self.class_1 = class_1
        self.class_2 = class_2
        self.k = k
        self.r0 = r0
        super().__init__()

    def gen_xml(self):
        return f'  <Bond class1="{self.class_1}" class2="{self.class_2}" length="{self.r0}" k="{self.k}"/>\n'


class HarmonicAngleForce(OpenMMXMLField):
    xml_field = "HarmonicAngleForce"

    def __init__(self, class_1, class_2, class_3, k, theta0):
        self.class_1 = class_1
        self.class_2 = class_2
        self.class_3 = class_3
        self.k = k
        self.theta0 = theta0
        super().__init__()

    def gen_xml(self):
        return f'  <Angle class1="{self.class_1}" class2="{self.class_2}" class3="{self.class_3}" angle="{self.theta0}" k="{self.k}"/>\n'


class PeriodicTorsionForceImproper(OpenMMXMLField):
    xml_field = "PeriodicTorsionForce"
    # HERE BE DRAGONS
    # Amber defines the central atom as the 3rd
    # OpenMM XML defines the central atom ast the 1st

    def __init__(self, class_1, class_2, class_3, class_4, k, phase, periodicity):
        self.class_1 = class_3  # Note the swap
        self.class_2 = class_2
        self.class_3 = class_1  # Note the swap
        self.class_4 = class_4
        self.k = k
        self.phase = phase
        self.periodicity = periodicity
        super().__init__()

    def gen_xml(self):
        return f'  <Improper type1="{self.class_1}" type2="{self.class_2}" type3="{self.class_3}" type4="{self.class_4}" k1="{self.k}" periodicity1="{int(self.periodicity)}" phase1="{self.phase}"/>\n'


class PeriodicTorsionForce(PeriodicTorsionForceImproper):
    def gen_xml(self):
        return f'  <Proper type1="{self.class_1}" type2="{self.class_2}" type3="{self.class_3}" type4="{self.class_4}" k1="{self.k}" periodicity1="{int(self.periodicity)}" phase1="{self.phase}"/>\n'


class NonbondedForce(OpenMMXMLField):
    xml_field = "NonbondedForce"

    def __init__(self, class_1, sigma, epsilon):
        self.class_1 = class_1
        self.sigma = sigma
        self.epsilon = epsilon
        self.charge = 0
        super().__init__()

    def gen_xml(self):
        return f'  <Atom type="{self.class_1}" charge="{self.charge}" sigma="{self.sigma}" epsilon="{self.epsilon}"/>\n'

    def gen_xml_header(self):
        #  TODO is this correct?
        return f' <{self.xml_field} coulomb14scale="0.833333" lj14scale="0.5">\n'


def get_section(key, line, f):
    if line == key["key"]:
        line = f.readline()  # Skip over header
        while line != "\n":
            parser = key["parser"]
            parser(line)
            line = f.readline()


def merge_phases(list_of_dihedrals):
    first = list_of_dihedrals[0]
    xml = f'  <Proper type1="{first.class_1}" type2="{first.class_2}" type3="{first.class_3}" type4="{first.class_4}"'
    for period, dihedral in enumerate(list_of_dihedrals, start=1):
        xml += f' k{period}="{dihedral.k}" periodicity{period}="{int(abs(dihedral.periodicity))}" phase{period}="{dihedral.phase}"'
    return xml + "/>\n"


def write_xml(xml_data, xml_name):
    with open(xml_name, "w") as f:
        f.write(xml_data)


if __name__ == "__main__":
    # ff_params = "eh-idtbr-frcmod"
    ff_params = "all-frcmod"
    xml = "<ForceField>\n"
    with open(ff_params) as f:
        line = f.readline()
        while line:
            for key, value in sections.items():
                get_section(value, line, f)
            line = f.readline()
    my_openMM = OpenMMXMLField()
    for field in my_openMM.open_xml_prams:
        field_items = my_openMM.open_xml_prams[field]
        idx = 0
        # Make header
        # Skip making the improper header
        if type(field_items[idx]) != PeriodicTorsionForceImproper:
            xml += field_items[0].gen_xml_header()
        while idx < len(field_items):
            if type(field_items[idx]) == PeriodicTorsionForce:
                if field_items[idx].periodicity < 0:
                    list_to_merge = []
                    while field_items[idx].periodicity < 0:
                        list_to_merge.append(deepcopy(field_items[idx]))
                        idx += 1
                    list_to_merge.append(deepcopy(field_items[idx]))
                    xml += merge_phases(list_to_merge)
                else:
                    xml += field_items[idx].gen_xml()
            else:
                xml += field_items[idx].gen_xml()
            idx += 1
        else:
            # skip making the period torsin footer
            if type(field_items[idx - 1]) != PeriodicTorsionForce:
                xml += field_items[idx - 1].gen_xml_closer()
    xml += "</ForceField>\n"
    write_xml(xml, "gaff.4fxml")
