# PlanckTon

## How to use

### Install

This will get less complicated when mBuild merges in a few PRs

```
git clone git@bitbucket.org:cmelab/planckton.git
cd planckton
conda create -n planckton
conda activate planckton
conda install -y --only-deps -c omnia -c cmelab -c mosdef mbuild foyer python=3.5
conda remove -y mbuild
conda install -y -c omnia -c conda-forge openmm=7.2.2 hoomd numpy=1.15.2 python=3.6
pip install -r requirements.txt
pip install .
pip install gsd signac
conda install numpy=1.15.2
conda install -c omnia -c mosdef mbuild
pytest # Run tests
```

You will also need to install [hoomd-blue](https://hoomd-blue.readthedocs.io/en/stable/)

### Run simulations

See example in `tests/test_sim.py`

Also see [planckton-flow](https://bitbucket.org/cmelab/planckton-flow)

### Add new compounds 


## How to develop

* Fork repo
* Add upstream `git remote add upstream git@bitbucket.org:cmelab/planckton.git`
* Code
* Submit PR

## Things to add

* Make the jobs easier to restart
** run up to
** don't re-init
** https://hoomd-blue.readthedocs.io/en/stable/restartable-jobs.html#temperature-ramp

## Debug notes 

`singularity exec --nv --bind $(pwd):/run/user/ planckton-test.simg python planckton/init.py`

`singularity exec --nv --bind $(pwd):/run/user/ planckton-test.simg python planckton/sim.py`

This repository is built to enable large sweeps exploring the self-assembly and charge transport 
for mixtures of organic photovoltaic compounds under various conditions.
To conduct these sweeps, multiple simulation tools are tied together 
to reproducibly and accurately generate these structures.

The simulations tools used here are:

* **mbuild** - builds compounds and initializes simulations

* **foyer** - selects the force-fields for the compounds

* **hoomd** - runs the simulations

* **signac** - manages and organizes the many variables that are run

* **signac-flow** - manages workflows

* **MorphCT** - determines charge transport through the morphology
