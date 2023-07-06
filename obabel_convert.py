from openbabel import openbabel

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("g09", "gjf")

mol = openbabel.OBMol()
obConversion.ReadFile(mol, "output11469.log")   # Open Babel will uncompress automatically

mol.AddHydrogens()

print(mol.NumAtoms())
print(mol.NumBonds())
print(mol.NumResidues())

obConversion.WriteFile(mol, '1abc.gjf')