from os import path
import warnings

import foyer.forcefields as ff
from foyer import Forcefield

FF_DIR = path.abspath(path.dirname(__file__))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    FORCE_FIELD = {
            "opv_gaff": Forcefield(path.join(FF_DIR, "gaff/opv_gaff.xml")),
            "oplsua-custom": Forcefield(
                path.join(FF_DIR, "oplsua/opls-custom.xml")
                ),
            "gaff": ff.load_GAFF(),
        }
