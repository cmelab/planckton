#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

AMBER_PREFIX="/home/mike/Projects/amber18"
PATH="${AMBER_PREFIX}/bin:${PATH}"

# Convert cml into mol2
INPUT_CMLS="$*"
for FILE in $INPUT_CMLS
do
    COMPOUND=${FILE%.cml}
    MOL="$COMPOUND".mol2
    MOL_TYPED="$COMPOUND"_typed.mol2
    FRCMOD="$COMPOUND"-frcmod
    echo Working on "$COMPOUND"
    echo Converting to mol2
    babel -i cml "$FILE"  -o mol2 "$MOL" # Convert cml into mol2
    echo cml to mol2 finished
    echo typing atoms
    antechamber -i "$MOL" -fi mol2 -o "$MOL_TYPED" -fo mol2 -c dc -pf y # Get atom types from amber
    echo atoms typed
    echo setting ff
    parmchk2 -i "$MOL_TYPED" -f mol2 -o "$FRCMOD" -a Y # Write out all the ff prams
    echo ff finished
done
