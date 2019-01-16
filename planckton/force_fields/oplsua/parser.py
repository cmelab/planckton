from dataclasses import dataclass, InitVar
from mass import MASS, NON_ELEMENT

PI = 3.141_592_653_589

FOYER_XML = True
# MIXED UA (still AA dihedrals) AA
FF_XML_TYPE = "MIXED"

if FF_XML_TYPE == "UA":
    ff_par = "oplsua.par.edits"
    ff_par_HEADER = 2
    ff_par_OPLS_TYPE_END = 434

    ff_sb = "oplsua.sb.edits"
    ff_sb_HEADER = 2
    ff_sb_HARMONIC_BOND_END = 178

    ff_sb_HARMONIC_ANGLE_START = 179

    ff_dihedrial = "../oplsaa/oplsaa.par.edits"
    ff_dihedrial_SKIP = 4038


elif FF_XML_TYPE == "AA":
    ff_par = "../oplsaa/oplsaa.par.edits"
    ff_par_HEADER = 2
    ff_par_OPLS_TYPE_END = 4031

    ff_sb = "../oplsaa/oplsaa.sb.edits"
    ff_sb_HEADER = 1
    ff_sb_HARMONIC_BOND_END = 451

    ff_sb_HARMONIC_ANGLE_START = 454

    ff_dihedrial = "../oplsaa/oplsaa.par.edits"
    ff_dihedrial_SKIP = 4038

elif FF_XML_TYPE == "MIXED":
    ff_par = "oplsua.par.edits"
    ff_par_HEADER = 2
    ff_par_OPLS_TYPE_END = 435

    ff_sb = "../oplsaa/oplsaa.sb.edits"
    ff_sb_HEADER = 1
    ff_sb_HARMONIC_BOND_END = 451

    ff_sb_HARMONIC_ANGLE_START = 454

    ff_dihedrial = "../oplsaa/oplsaa.par.edits"
    ff_dihedrial_SKIP = 4038


# Since we have a few PRs to make for this to work in Foyer I don't want to do
# Anything by hand, so I will make a dictionary to fill out params we need

OPLS_INFO = {
    "opls_068": {
        "def": "[C!R;X1][C;!R]",
        "desc": "CH3 (C2) N-ALKANES",
        "doi": "10.1021/ja00334a030",
    },
    "opls_071": {
        "def": "[C!R;X2][C;!R]C",
        "desc": "CH2 (SP3) ALKANES",
        "doi": "10.1021/ja00334a030",
    },
    "opls_230": {
        "def": "[CR;X2][C;R][C;R]",
        "desc": "Benzene C - 12 site JACS,112,4768-90",
        "doi": "10.1021/ja00168a022",
    },
    "opls_380": {"def": "[CR;X2][C;R][S;R]", "desc": "CE1 in HID,HIE", "doi": ""},
    "opls_500": {
        "def": "[SR;X2][C;R][S;R]",
        "desc": "Taken from AA opls_2003",
        "doi": "",
    },
}


@dataclass(frozen=False)
class OPLSUA_type:
    opls_type: str
    atomic_number: int
    atomic_name: str
    atomic_class: str
    charge: float
    sigma: float
    epsilon: float
    mass_dic: InitVar
    info_dic: InitVar
    element_translate: InitVar
    _def = ""
    _desc = ""
    _doi = ""

    def __post_init__(self, mass_dic, info_dic, element_translate):
        self.mass = mass_dic[self.atomic_name]  # Need to make a mass dic for AA
        try:
            info = info_dic[self.opls_type]
            self._def = info["def"]
            self._desc = info["desc"]
            self._doi = info["doi"]
        except KeyError:
            pass
        try:
            self.atomic_name = element_translate[self.atomic_name]
        except KeyError:
            self.atomic_name = self.atomic_name
            pass


@dataclass(frozen=True)
class HarmonicBond:
    class_1: str
    class_2: str
    k: float
    lenght: float


@dataclass(frozen=True)
class HarmonicAngle:
    class_1: str
    class_2: str
    class_3: str
    k: float
    angle: float


@dataclass(frozen=True)
class RBTorsion:
    class_1: str
    class_2: str
    class_3: str
    class_4: str
    c0: float
    c1: float
    c2: float
    c3: float
    c4: float
    c5: float


def write_xml(xml_data, xml_name):
    with open(xml_name, "w") as f:
        f.write(xml_data)


def opls_dihedral_to_RB_torsion(f1, f2, f3, f4):
    c0 = f2 + 0.5 * (f1 + f3)
    c1 = 0.5 * (-f1 + 3 * f3)
    c2 = -f2 + 4 * f4
    c3 = -2 * f3
    c4 = -4 * f4
    c5 = 0
    return c0, c1, c2, c3, c4, c5


oplsua_list = []
problem_lines = []


with open(ff_par, "r") as f:
    ff_par_lines = f.readlines()

for line in ff_par_lines[ff_par_HEADER:ff_par_OPLS_TYPE_END]:
    raw_opls_type = line.strip(" ").strip().split(" ")
    opls_param_array = list(filter(None, raw_opls_type))
    if opls_param_array[0][0] == "#":
        problem_lines.append(opls_param_array)
        continue
    opls_type = f"opls_{int(opls_param_array[0]):03d}"
    try:
        atomic_name = opls_param_array[2]
        atomic_class = atomic_name
    except IndexError:
        problem_lines.append(opls_param_array)
        continue
    # Skip the ?? and DU types
    if any(_ in atomic_name for _ in ("??", "DU")):
        problem_lines.append(opls_param_array)
        continue
    atomic_number = opls_param_array[1]
    charge = opls_param_array[3]
    # We need to convert Å to nm
    sigma = float(opls_param_array[4]) / 10
    # We need to convert kcal/mol to kJ/mol
    epsilon = float(opls_param_array[5]) * 4.184
    new_opls_type = OPLSUA_type(
        opls_type,
        atomic_number,
        atomic_name,
        atomic_class,
        charge,
        sigma,
        epsilon,
        mass_dic=MASS,
        info_dic=OPLS_INFO,
        element_translate=NON_ELEMENT,
    )
    oplsua_list.append(new_opls_type)

print(f"lines skipped in {ff_par} {problem_lines} \n")

with open(ff_sb, "r") as f:
    ff_sb_lines = f.readlines()

problem_lines = []
harmonic_bond_types = []
for line in ff_sb_lines[ff_sb_HEADER:ff_sb_HARMONIC_BOND_END]:
    raw_harmonic_bond_type = line.strip(" ").strip().split(" ")
    if line.startswith(("#", "*", '"')):
        problem_lines.append(line)
        continue
    harmonic_bond_array = list(filter(None, raw_harmonic_bond_type))
    print(harmonic_bond_array)
    # In some cases, there is a space between two classes, ie C -CA this causes an
    # issue when we split things up, so we check to see if we need to concatenate
    # the first two items in our list
    if "-" not in harmonic_bond_array[0]:
        harmonic_bond_array[:2] = ["".join(harmonic_bond_array[:2])]
    class_1, class_2 = harmonic_bond_array[0].split("-")
    # We need to convert kcal/(Å**2 mol) to kJ/(nm**2 mol)
    k = float(harmonic_bond_array[1]) * (4.184 * 100 * 2)
    # We need to convert Å to nm
    lenght = float(harmonic_bond_array[2]) / 10
    new_harmonic_bond = HarmonicBond(class_1, class_2, k, lenght)
    harmonic_bond_types.append(new_harmonic_bond)

harmonic_angle_types = []
for line in ff_sb_lines[ff_sb_HARMONIC_ANGLE_START:]:
    # ** seems to indicate a comment
    if line.startswith("**"):
        problem_lines.append(line)
        continue
    raw_harmonic_angle_type = line.strip(" ").strip().split(" ")
    harmonic_angle_array = list(filter(None, raw_harmonic_angle_type))
    while len(harmonic_angle_array[0].split("-")) != 3:
        harmonic_angle_array[0] += harmonic_angle_array.pop(1)
    class_1, class_2, class_3 = harmonic_angle_array[0].split("-")
    # We need to convert kcal/(deg**2 mol) to kJ/(rad**2 mol)
    k = float(harmonic_angle_array[1]) * (4.184 * 2)
    # We need to converd deg to rad
    angle = float(harmonic_angle_array[2]) * (PI / 180)
    new_harmonic_angle = HarmonicAngle(class_1, class_2, class_3, k, angle)
    harmonic_angle_types.append(new_harmonic_angle)

print(f"lines skipped in {ff_dihedrial} {problem_lines} \n")

with open(ff_dihedrial, "r") as f:
    ff_dihedrial_lines = f.readlines()
rb_torsion_list = []
problem_lines = []

for line in ff_dihedrial_lines[ff_dihedrial_SKIP:]:
    raw_opls_dihedrial_type = line.strip(" ").strip().split(" ")
    dihedrial_array = list(filter(None, raw_opls_dihedrial_type))
    # Check to see if line is a dummy line or leave blank
    if any(_ in dihedrial_array for _ in ("Dummy", "UNASSIGNED")):
        problem_lines.append(dihedrial_array)
        continue
    # Like with the angles and bonds, we need to deal with spaces
    # But with dihedrials, the index is different
    # We also need to catch bad lines
    try:
        while len(dihedrial_array[5].split("-")) != 4:
            dihedrial_array[5] += dihedrial_array.pop(6)
        # We need to remove the double bond notation "="
        class_1, class_2, class_3, class_4 = (
            dihedrial_array[5].replace("=", "").split("-")
        )
    except IndexError:
        problem_lines.append(dihedrial_array)
        continue
    f1, f2, f3, f4 = map(float, dihedrial_array[1:5])
    # Convert to RB style torsions
    c_coefs = opls_dihedral_to_RB_torsion(f1, f2, f3, f4)
    new_rb_torsion = RBTorsion(class_1, class_2, class_3, class_4, *c_coefs)
    rb_torsion_list.append(new_rb_torsion)

print(f"lines skipped in {ff_dihedrial} {problem_lines} \n")


openMM_xml = "<ForceField>\n"

# Atom Types
openMM_xml += "<AtomTypes>\n"
for opls_type in oplsua_list:
    if FOYER_XML:
        # Check if def is empty (https://github.com/mosdef-hub/foyer/issues/179)
        if not opls_type._def:
            openMM_xml += f' <Type name="{opls_type.opls_type}" class="{opls_type.atomic_class}" element="{opls_type.atomic_name}" mass="{opls_type.mass}"/>\n'
        else:
            openMM_xml += f' <Type name="{opls_type.opls_type}" class="{opls_type.atomic_class}" element="{opls_type.atomic_name}" mass="{opls_type.mass}" def="{opls_type._def}" desc="{opls_type._desc}" doi="{opls_type._doi}"/>\n'
    else:
        openMM_xml += f' <Type name="{opls_type.opls_type}" class="{opls_type.atomic_class}" element="{opls_type.atomic_name}" mass="{opls_type.mass}"/>\n'
openMM_xml += "</AtomTypes>\n"

# Harmonic Bond Force
openMM_xml += "<HarmonicBondForce>\n"
for harmonic_bond in harmonic_bond_types:
    openMM_xml += f' <Bond class1="{harmonic_bond.class_1}" class2="{harmonic_bond.class_2}" length="{harmonic_bond.lenght}" k="{harmonic_bond.k}"/>\n'
openMM_xml += "</HarmonicBondForce>\n"

# Harmonic Angle Force
openMM_xml += "<HarmonicAngleForce>\n"
for harmonic_angle in harmonic_angle_types:
    openMM_xml += f' <Angle class1="{harmonic_angle.class_1}" class2="{harmonic_angle.class_2}" class3="{harmonic_angle.class_3}" angle="{harmonic_angle.angle}" k="{harmonic_angle.k}"/>\n'
openMM_xml += "</HarmonicAngleForce>\n"

# RB Torsionorsion Force
openMM_xml += "<RBTorsionForce>\n"
for rb_torsion in rb_torsion_list:
    openMM_xml += f' <Proper class1="{rb_torsion.class_1}" class2="{rb_torsion.class_2}" class3="{rb_torsion.class_3}" class4="{rb_torsion.class_4}" c0="{rb_torsion.c0}" c1="{rb_torsion.c1}" c2="{rb_torsion.c2}" c3="{rb_torsion.c3}" c4="{rb_torsion.c4}" c5="{rb_torsion.c5}"/>\n'
openMM_xml += "</RBTorsionForce>\n"
# Non-bonded Force
openMM_xml += '<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">\n'
for opls_type in oplsua_list:
    openMM_xml += f' <Atom type="{opls_type.opls_type}" charge="{opls_type.charge}" sigma="{opls_type.sigma}" epsilon="{opls_type.epsilon}"/>\n'
openMM_xml += "</NonbondedForce>\n"

openMM_xml += "</ForceField>\n"
write_xml(openMM_xml, f"opls-{FF_XML_TYPE}.xml")
