# This repository is no longer under active development and has been archived.
See [**flowerMD**](https://github.com/cmelab/flowerMD) for an updated package that performs similar workflows.


PlanckTon enables exploration of the self-assembly of organic photovoltaic (OPV) compound mixtures under various conditions.
Multiple simulation tools are tied together to reproducibly and accurately generate these structures.
For managing large parameter spaces and submitting jobs to clusters, we recommend [PlanckTon-flow](https://github.com/cmelab/planckton-flow) which leverages the [Signac](https://docs.signac.io/en/latest/) framework.

The simulation tools used by PlanckTon are:

* [**mBuild**](https://github.com/mosdef-hub/mbuild) - builds compounds and initializes simulations

* [**foyer**](https://foyer.mosdef.org/en/stable/) - manages force-field information

* [**HOOMD-blue**](https://hoomd-blue.readthedocs.io/en/latest/) - runs the simulations
