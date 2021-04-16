# OPLS-UA Parsing

## Conversion methods
* To convert angstroms to nano meters, a conversion factor of 1/10 is used.
* To convert kcal/mol to kJ/mol, a conversion factor of 4.184 is used.
* To convert kcal/(Å**2 mol) to kJ/(nm**2 mol), a conversion factor of 4.184 * 100 is used.
* To convert radians to degrees, a conversion factor of 3.141592653589/180 is used.
* OPLS style dihedrals are converted to RB torsions with the following method:
  ```
  c0 = f2+0.5*(f1+f3)
  c1 = 0.5*(-f1+3*f3)
  c2 = -f2 + 4*f4
  c3 = -2*f3
  c4 = -4*f4
  c5 = 0
  ```
  They can be converted back to OPLS style with the following method:
  ```
  f1 = (-1.5 * c3) - (2 * c1)
  f2 = c0 + c1 + c3
  f3 = -0.5 * c3
  f4 = -0.25 * c4
  ```

## Notes
* Line 155 in oplsua.sb has a duplicate param, but different K and r, not sure how to handle this, in oplsua.sb.edits, I've renamed the HO to HM to reflect the authors last name, but nothing will ever be typed HM.
* `oplsua.xml` Has UA types, bonds, and angles but AA dihedrals from `oplsaa.par.edits`
* Some atom types in `oplsua.par` are `??` The parsing code skips these types.
* Dihedirals can have wild cards by specifying and empty class/type, but bonds and angles cannot. Will need to check to see if any wild cards in are in the xml.
* Something to think about with the mass dictionary, are atom classes reused? I don't know.
* Hmmm I do think there is reuse, the net effect is I'll be off by a mass of H (or two).
* The more I have read about OPLS UA and AA, the more I think we shouldn't use UA for our actual sweep. We should look at amber or other CG FF, it looks like Jorgensen has moved away from the UA
* Maybe as an aside, we do some AA sims and paramterize a new OPV CG FF? Skip UA and go onto MS-IBI?
* I need to make the AA stuff parsed into an openMM xml as well so I can grab AA params more easily 
* I want to test using this site for AA http://erg.biophys.msu.ru/erg/tpp/ then droping H and see what happens
* I need to check units on dihedrals!
