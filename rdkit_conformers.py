## Conformer  generating script for the Tresca Lab
## Version 1.0 by Blakely Tresca, 2023
## Original Code from RDKit Cookbook and mcsorkun

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# read spreadsheet title from command lines
#req_list = sys.argv[1]

# ask for the input structure
inp_str = input("Input structure as mol file or SMILES: ")
print(inp_str)
if inp_str == 'n':
    exit(0)
elif inp_str.find('.mol') != -1:
    m = Chem.MolFromMolFile(inp_str)
    m.SetProp("_Name",inp_str.replace('.mol', ''))
else:
    m = Chem.MolFromSmiles(inp_str)

if m is None:
    print('Import fail, check the .mol file')
    exit(0)

#print("Hello World")
print(Chem.MolToSmiles(m))
AllChem.Compute2DCoords(m)
print(Chem.MolToMolBlock(m))

#add hydrogens
m_Hs = Chem.AddHs(m)

#generate a 3D structure
m_3d = m_Hs
AllChem.EmbedMolecule(m_3d,randomSeed=0xf00d)
print(Chem.MolToMolBlock(m_3d),file=open(inp_str.replace('.mol', '_3d.mol'),'w+'))

#save an image of the 2D structure
Draw.MolToFile(m,'molecule.png')

#generate a conformer library
rmslist = []
#params = AllChem.ETKDGv2()
m_confs = AllChem.EmbedMultipleConfs(m_3d, numConfs=100, numThreads=0)
AllChem.AlignMolConformers(m_3d, RMSlist=rmslist)
rms = AllChem.GetConformerRMS(m_3d, 1, 9, prealigned=True)
res = AllChem.MMFFOptimizeMoleculeConfs(m_3d, numThreads=0)

open(inp_str.replace('.mol', '_conf.sdf'),'w+')

#save the library to individual mol files
for m_conf in m_confs:
    print(m_conf)
    #print(Chem.MolToMolBlock(m_3d, confId=m_conf))
    #print(Chem.MolToMolBlock(m_3d, confId=m_conf),file=open('conf_'+str(m_conf)+'.mol','a+'))
    # (needs troubleshooting) m_3d.SetProp("_Name",inp_str.replace('.mol', f"conf {m_conf}"), confId=m_conf)
    print(f'conf_{m_conf} {Chem.MolToMolBlock(m_3d, confId=m_conf)}',file=open(inp_str.replace('.mol', '_conf.sdf'),'a+'))
    print(f'\n{res[m_conf]}\n$$$$',file=open(inp_str.replace('.mol', '_conf.sdf'),'a+'))
    #Draw.MolToFile(m_3d,'molecule.png')

#mols = [m_conf for m_conf in m_3d]
#img=Draw.MolsToGridImage(m_3d[0:4],molsPerRow=4,subImgSize=(200,200), legends=None)
#img.save('molgrid.png')
