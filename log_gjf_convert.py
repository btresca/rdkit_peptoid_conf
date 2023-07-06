# https://github.com/rdkit/rdkit/issues/3310

from openbabel import openbabel
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def add_nitrogen_charges(bad_mol):
    bad_mol.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(bad_mol)
    if not ps:
        Chem.SanitizeMol(bad_mol)
        return bad_mol
    for p in ps:
        if p.GetType()=='AtomValenceException':
            at = bad_mol.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                at.SetFormalCharge(1)
    Chem.SanitizeMol(bad_mol)
    return bad_mol

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("g09", "mol")

#mol = openbabel.OBMol()
#obConversion.ReadFile(mol, "output11469.log")   # Open Babel will uncompress automatically

#print(mol.NumAtoms())
#print(mol.NumBonds())
#print(mol.NumResidues())
#print(mol.GetTotalCharge())

for molecule in pybel.readfile("g09","output11469.log"):
    print(molecule.molwt)
    print(molecule.charge)
    mol = molecule

# obConversion.WriteFile(mol, '1abc.gjf')

#outMDL = obConversion.WriteString(mol)
#print(outMDL)

outMDL = mol.write("mol")
print(outMDL)

#ms = Chem.SDMolSupplier()
#ms.SetData(outMDL)
#rd_mol = next(ms)
rd_mol = Chem.MolFromMolBlock(outMDL,sanitize=False)
rd_mol = add_nitrogen_charges(rd_mol)
print(rd_mol)
Draw.MolToFile(rd_mol,'molecule.png')
print(Chem.MolToMolBlock(rd_mol),file=open('molecule.mol','w+'))




print("Done")