import os
import glob
from pathlib import Path

COMPOUND_DIR = Path(__file__).parent / "compounds"
COMPOUND_FILE = {}
for compound_file in glob.glob(os.path.join(COMPOUND_DIR, "*_typed.mol2")):
    # The the last part of the path, drop the _typed.mol2
    compound_name = compound_file.split("/")[-1].split("_typed.mol2")[0]
    COMPOUND_FILE[compound_name] = compound_file
