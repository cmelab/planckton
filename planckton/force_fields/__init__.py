import os

from foyer import Forcefield

FF_DIR = os.path.abspath(os.path.dirname(__file__))
FORCE_FIELD = {
        "opv_gaff": Forcefield(
            os.path.join(FF_DIR, "gaff/opv_gaff.xml")
            ),
        "oplsua-custom": Forcefield(
            os.path.join(FF_DIR, "oplsua/opls-custom.xml")
            ),
        }
