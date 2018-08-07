# PlanckTon

This repository is built to enable large sweeps exploring the self-assembly and charge transport 
for mixtures of organic photovoltaic compounds under various conditions.
To conduct these sweeps, multiple simulation tools are tied together 
to reproducably and accurately generate these structures.

The simulations tools used here are:

* **mbuild** - builds compounds and initializes simulations

* **foyer** - selects the force-fields for the compounds

* **hoomd** - runs the simulations

* **signac** - manages and organizes the many variables that are run

* **signac-flow** - manages workflows

* **MorphCT** - determines charge tranport through the morphology