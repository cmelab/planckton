# PlanckTon
[![test](https://github.com/cmelab/planckton/workflows/test/badge.svg)](https://github.com/cmelab/planckton/actions?query=workflow%3Atest)
[![build](https://github.com/cmelab/planckton/workflows/build/badge.svg)](https://github.com/cmelab/planckton/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/cmelab/planckton/branch/master/graph/badge.svg?token=5KYVHWMT28)](https://codecov.io/gh/cmelab/planckton/)
[![License](https://img.shields.io/badge/license-GPLv3-green.svg)](LICENSE.md)
[![Contributors](https://img.shields.io/github/contributors-anon/cmelab/planckton.svg?style=flat)](https://github.com/cmelab/planckton/graphs/contributors)

PlanckTon enables exploration of the self-assembly of organic photovoltaic compound mixtures under various conditions.
Multiple simulation tools are tied together to reproducibly and accurately generate these structures.

The simulation tools used by PlanckTon are:

* [**mBuild**](https://github.com/mosdef-hub/mbuild) - builds compounds and initializes simulations

* [**foyer**](https://foyer.mosdef.org/en/stable/) - manages force-field information

* [**HOOMD-blue**](https://hoomd-blue.readthedocs.io/en/latest/) - runs the simulations


## How to use

### Install
#### Using a container
To use PlanckTon in a prebuilt container (using [Singularity](https://singularity.lbl.gov/)), run:
```bash
singularity pull docker://cmelab/planckton_cpu:latest
singularity exec planckton_cpu_latest.sif bash
```

**Or** using [Docker](https://docs.docker.com/), run:
```bash
docker pull cmelab/planckton_cpu:latest
docker run -it cmelab/planckton_cpu:latest
```

#### Custom install
To create a local environment with [conda](https://docs.conda.io/en/latest/miniconda.html), run:
```bash
conda env create -f environment.yml
conda activate planckton
```
And to test your installation, run:
```
pytest
```

### Run simulations

See example in `tests/test_sim.py`

Also see [planckton-flow](https://github.com/cmelab/planckton-flow)

## How to develop

* Fork and clone this repo ([how to make a fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo))
* Add this repository as upstream `git remote add upstream git@github.com:cmelab/planckton.git`
* Modify the code and push to your fork
* Submit a pull request (PR) ([how to create a PR](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork))
