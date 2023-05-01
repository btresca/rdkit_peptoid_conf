# rdkit_peptoid_conf
Python scripts to generate conformers using RDKit. \n
Install RDKit in python environment following Conda instructions here https://www.rdkit.org/docs/Install.html \n
Input can be .mol file or SMILES, either can be generated from ChemDoodle. \n
Output files include: \n
  image.png - 2D image of molecule
  mol_3D.mol - initial 3D coordinates of molecule with H's
  mol_min-E-conf-MMFF.sdf - 3D coordinates of minimum Ener conformer by MMFF optimization
  mol_conf.sdf - 3D coordiantes of all conformers meeting criteria
  mol_E_list.csv - list of conformer Index and rel. E
