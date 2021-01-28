import os
import glob

from foyer import Forcefield

FORCE_FIELD_DIR = os.path.abspath(os.path.dirname(__file__))
FORCE_FIELD = {}
for force_field_file in glob.glob(os.path.join(FORCE_FIELD_DIR, "*/*.xml")):
    force_field_name = force_field_file.split("/")[-1].split(".xml")[0]
    FORCE_FIELD[force_field_name] = Forcefield(force_field_file)
