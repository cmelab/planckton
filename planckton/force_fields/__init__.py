from os.path import join
import warnings

from foyer import Forcefield

FF_DIR = os.path.abspath(os.path.dirname(__file__))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    FORCE_FIELD = {
            "opv_gaff": Forcefield(join(FF_DIR, "gaff/opv_gaff.xml")),
            "oplsua-custom": Forcefield(join(FF_DIR, "oplsua/opls-custom.xml")),
        }
