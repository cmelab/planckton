import os
import glob

COMPOUND_DIR = os.path.abspath(os.path.dirname(__file__))
COMPOUND_FILE = {}
for compound_file in glob.glob(os.path.join(COMPOUND_DIR, "*_typed.mol2")):
    # The the last part of the path, drop the _typed.mol2
    compound_name = compound_file.split("/")[-1].split("_typed.mol2")[0]
    COMPOUND_FILE[compound_name] = compound_file
