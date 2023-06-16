## Conformer generating script for the Tresca Lab
## Version 2.0 by Blakely Tresca, 2023
## Original code from RDKit Cookbook and mcsorkun https://github.com/mcsorkun/Conformer-Search

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import math
import os

# read argument from command line
#req_list = sys.argv[1]

# Number of conformers to be generated
num_of_conformer = 500
max_iter = 500
rms_thresh = 0.5
# Default values for min energy conformer
min_energy_UFF = 10000
min_energy_index_UFF = 0
min_energy_MMFF = 10000
min_energy_index_MMFF = 0
# Default values for keeping structures
max_ener = 40  # kcal/mol
min_conf = 20
conf_keep = []
conf_ener = []
conf_list = 0

# Define method to find dihedral angles
def mes_dihedral():
    # Prompt the user to define the atom indices for dihedral calculation
    atom1 = int(input("Enter the atom index for Calpha_2: "))
    atom2 = int(input("Enter the atom index for C=O_2: "))
    atom3 = int(input("Enter the atom index for N_1: "))
    atom4 = int(input("Enter the atom index for Calpha_1: "))
    atom5 = int(input("Enter the atom index for C=O_1: "))
    atom6 = int(input("Enter the atom index for N_2: "))
    if inp_length >= 2:
        atom7 = int(input("Enter the atom index for N_2: "))
        atom8 = int(input("Enter the atom index for N_2: "))

# ask for the input structure
inp_str = input("Input structure as mol file or SMILES: ")
print(inp_str)
if inp_str == 'n':
    exit(0)
elif inp_str.find('.mol') != -1:
    m = Chem.MolFromMolFile(inp_str)
    m.SetProp("_Name", inp_str.replace('.mol', ''))
    mol_name = inp_str.replace('.mol', '')
else:
    m = Chem.MolFromSmiles(inp_str)
    mol_name = input("What shall we call this molecule: ")
    print(mol_name)
    m.SetProp("_Name", mol_name)

if m is None:
    print('Import fail, check the .mol file')
    exit(0)

inp_dihedral = input("Measure dihedrals? (y/n) ")
print(inp_str)

if inp_dihedral == 'y':
    inp_length = int(input("How many amides with dihedrals? "))
    print(inp_length)

#print("Hello World")
print(Chem.MolToSmiles(m))
AllChem.Compute2DCoords(m)
print(Chem.MolToMolBlock(m))

#add hydrogens
m_Hs = Chem.AddHs(m)

# Make directory for output files
if not os.path.exists(mol_name):
    os.mkdir(mol_name)
output_dir = f'{os.path.join(os.getcwd(), mol_name)}'
# Create the output folder
# output_folder = os.path.join(os.getcwd(), mol_name)
# os.makedirs(output_folder, exist_ok=True)



#generate a 3D structure
m_3d = m_Hs
AllChem.EmbedMolecule(m_3d, randomSeed=0xf00d)
print(Chem.MolToMolBlock(m_3d), file=open(os.path.join(output_dir, mol_name + '_3d.mol'), 'w+'))

for atom in m.GetAtoms():
    # For each atom, set the property "atomNote" to a index+1 of the atom
    atom.SetProp("atomNote", str(atom.GetIdx()+1))

#save an image of the 2D structure
Draw.MolToFile(m, os.path.join(output_dir, mol_name + '.png'), size=(800, 800))

#generate a conformer library
# Generate conformers (stored inside the mol object)
print("\nFeeding hamsters")
rmslist = []
m_confs = AllChem.EmbedMultipleConfs(m_3d, numConfs=num_of_conformer, pruneRmsThresh=rms_thresh, numThreads=0)
conf_ids = list(m_confs)  # create list of conformer IDs

# Align molecules and optimize energy with MMFF
print("\nOptimizing conformers with MMFF")
AllChem.AlignMolConformers(m_3d, RMSlist=rmslist)
rms = AllChem.GetConformerRMS(m_3d, 1, 9, prealigned=True)
results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(m_3d, maxIters=max_iter, numThreads=0)

# Search for the min energy conformer from results(tuple(is_converged,energy))
print("\nSearching for min E conformer by MMFF ")   
for index, result in enumerate(results_MMFF):
    if(min_energy_MMFF > result[1]):
        min_energy_MMFF = result[1]
        min_energy_index_MMFF = index
        print(min_energy_index_MMFF, ":", min_energy_MMFF)

## deltaE = results_MMFF - min_energy_MMFF

# Write minimum energy conformer into a SDF file
w = Chem.SDWriter(os.path.join(output_dir, mol_name + '_min-E-conf-MMFF.sdf'))
w.write(Chem.Mol(m_3d,False,min_energy_index_MMFF))
w.flush()
w.close()

# Search for conformers within max_ener of min energy conformer
print("\nSearching for conformers to keep ")
for index, result in enumerate(results_MMFF):
    rel_ener = result[1] - min_energy_MMFF
    if (rel_ener < max_ener):
        # Remove duplicate conformers with deltaE less than abs_tol
        close_e = min(conf_ener, key = lambda x:abs(x-result[1]), default = 0)
        rel_close_e = close_e - min_energy_MMFF
        if not math.isclose(rel_ener, rel_close_e, abs_tol = 0.1):
            conf_keep.append(index)
            conf_ener.append(result[1])
            print(f'{index}:{result[1]}')

# Sort the conformers by energy (lowest to highest)
sorted_conf = sorted(zip(conf_keep, conf_ener), key=lambda x: x[1])
conf_keep, conf_ener = zip(*sorted_conf)

output_sdf = os.path.join(output_dir, mol_name + '_conf.sdf')
open(output_sdf, 'w+')
print("index, energy",file=open(os.path.join(output_dir, mol_name + '_E_list.csv'),'w+'))

for atom in m_3d.GetAtoms():
    # For each atom, set the property "atomNote" to a index+1 of the atom
    atom.SetProp("atomNote", str(atom.GetIdx()+1))

#save an image of the 3D structure
Draw.MolToFile(m_3d, os.path.join(output_dir, mol_name + '_3D.png'), size=(800, 800))

#save the library to individual mol files
for index in conf_keep:
    print(f'conf_{index} {Chem.MolToMolBlock(m_3d, confId=index)}', file=open(output_sdf,'a+'))
    rel_ener = conf_ener[conf_list] - min_energy_MMFF
    print(f'> <Relative Energy>\n{rel_ener}\n\n$$$$', file=open(output_sdf, 'a+'))
    print(f'{index}, {rel_ener}', file=open(os.path.join(output_dir, mol_name + '_E_list.csv'),'a+'))

    conf_list = conf_list + 1

#mols = [m_conf for m_conf in m_3d]
#img=Draw.MolsToGridImage(m_3d[0:4],molsPerRow=4,subImgSize=(200,200), legends=None)
#img.save('molgrid.png')
