"""Forcefields included with PlanckTon."""

import warnings
from os import path

import foyer
from foyer import Forcefield

FF_DIR = path.abspath(path.dirname(__file__))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    FORCEFIELD = {
        "gaff-custom": Forcefield(path.join(FF_DIR, "gaff/opv_gaff.xml")),
        "oplsua-custom": Forcefield(path.join(FF_DIR, "oplsua/oplsua.xml")),
        "gaff": foyer.forcefields.load_GAFF(),
    }
