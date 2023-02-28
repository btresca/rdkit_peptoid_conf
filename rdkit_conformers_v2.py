## Conformer  generating script for the Tresca Lab
## Version 2s.0 by Blakely Tresca, 2023
## Original Code from RDKit Cookbook and mcsorkun

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# read spreadsheet title from command lines
#req_list = sys.argv[1]

# Number of conformers to be generated
num_of_conformer=100
max_iter=500
# Default values for min energy conformer
min_energy_UFF=10000
min_energy_index_UFF=0
min_energy_MMFF=10000
min_energy_index_MMFF=0
# Default values for keeping structures
max_ener=40 #kcal/mol
min_conf=20

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
# Generate conformers (stored in side the mol object)
rmslist = []
params = AllChem.ETKDGv2()
m_confs = AllChem.EmbedMultipleConfs(m_3d, numConfs=1000, numThreads=0, params = params)
conf_ids = list(m_confs) #create list of conformer IDs

AllChem.AlignMolConformers(m_3d, RMSlist=rmslist)
rms = AllChem.GetConformerRMS(m_3d, 1, 9, prealigned=True)
results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(m_3d, maxIters=max_iter, numThreads=0)

# Search for the min energy conformer from results(tuple(is_converged,energy))
print("\nSearching conformers by MMFF ")   
for index, result in enumerate(results_MMFF):
    if(min_energy_MMFF>result[1]):       
        min_energy_MMFF=result[1]
        min_energy_index_MMFF=index
        print(min_energy_index_MMFF,":",min_energy_MMFF)

## deltaE = results_MMFF - min_energy_MMFF

# Write minimum energy conformer into a SDF file
w = Chem.SDWriter('minimum-energy-conformer-MMFF.sdf')
w.write(Chem.Mol(m_3d,False,min_energy_index_MMFF))
w.flush()  
w.close()

## To Do: Calc deltaE, create list of IDs:deltaE, pick lowest by parameters, write keep to sdf

open(inp_str.replace('.mol', '_conf.sdf'),'w+')

#save the library to individual mol files
for index, result in enumerate(results_MMFF):
    print(index,":",result)
    #print(Chem.MolToMolBlock(m_3d, confId=m_conf))
    #print(Chem.MolToMolBlock(m_3d, confId=m_conf),file=open('conf_'+str(m_conf)+'.mol','a+'))
    # (needs troubleshooting) m_3d.SetProp("_Name",inp_str.replace('.mol', f"conf {m_conf}"), confId=m_conf)
    print(Chem.MolToMolBlock(m_3d, confId=m_conf),file=open(inp_str.replace('.mol', '_conf.sdf'),'a+'))
    print('$$$$',file=open(inp_str.replace('.mol', '_conf.sdf'),'a+'))
    #Draw.MolToFile(m_3d,'molecule.png')

#mols = [m_conf for m_conf in m_3d]
#img=Draw.MolsToGridImage(m_3d[0:4],molsPerRow=4,subImgSize=(200,200), legends=None)
#img.save('molgrid.png')
