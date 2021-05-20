# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## 0.4.0 - 2021-04-20
- Code has been refactored:
    - all typed compounds have gaff added in their name so it is more clear what FF they are typed with (e.g. P3HT_typed.mol2 --> P3HT-gaff_typed.mol2)
    - force_fields directory renamed to forcefields
    - COMPOUND_FILE and FORCE_FIELD dicts renamed to COMPOUND and FORCEFIELD
    - deleted cml and untyped mol2 files in compounds directory
- Added pre-commit yaml to automatically lint the code:
    - many whitespace changes were made (eol and eof whitespace deleted)
    - docstrings updated

## 0.3.1 - 2021-04-16
- add error checking for unimplented functionality

## 0.3.0 - 2021-04-16
- begin supporting smiles/outside forcefields

## 0.2.1 - 2021-03-15
- restart.gsd tracks momenta
- velocities are not randomized on restart

## 0.2.0 - 2021-03-08
- add restart capability
- gsd files track image

## 0.1.5 - 2021-01-29
- use unyt and ele packages for handling units and elements
- use gsd instead of deprecated xml files
- use `create_hoomd_simulation` function from mbuild
- remove `cme_utils` dependency

## 0.0.2 - 2019-01-16
### Added
- starting compounds
- forcefield
- init code
- sim code
- lots more :)

## 0.0.1 - 2018-10-16
### Added
- Everything :)
