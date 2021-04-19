# How to add new compounds

## Set up

First we need a few more tools:

```
conda install -c openbabel -c ambermd ambertools=19 openbabel numpy
```

## Overview

The basic process is this:

1. Convert `cml` to `mol2` using babel
1. Pass `mol2` to `antechamber` and have it generate a `mol2` that is typed
1. Use `parmchk2` to write out the needed force field params baised on the atom typeing, this creates a `-frcmod` file.
1. `cat` all of the `-frcmod` files into a `all-frcmod` file, ie `cat *-frcmod >> all-frcmod`
1. `python parser.py` will then read in the `all-frcmod` file and generate a `gaff.4fxml` xml file that will work with foyer.
1. Then I check to see if we added any new atom types by `diff gaff.4fxml ../force_fields/gaff/opv_gaff.xml` and see if anything new is added.
1. If there is nothing new, then all we have to do is `git add COMPOUND_typed.mol2` and then we can access the compound in planckton.
1. If there is something new, `cp gaff.4fxml ../force_fields/gaff/opv_gaff.xml` and add the compound like we did above, and commit the results.
