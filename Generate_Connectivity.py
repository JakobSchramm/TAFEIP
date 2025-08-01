import numpy as np
from ase.io import read
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs, get_connectivity_matrix
from rdkit import Chem
import ast

def con_mat_from_struc(name): # connectivity matrix of carbon structure from xyz file
    mol = read(name)
    mol_nl = NeighborList(natural_cutoffs(mol), self_interaction=False, bothways=True)
    mol_nl.update(mol)
    carbon = Atoms()
    for ind,atom in enumerate(mol):
        if atom.symbol == "C":
            if len(mol_nl.get_neighbors(ind)[0]) != 3:
                raise SystemExit("Error, one carbon is not sp2 hybridized!")
            carbon.append(atom)
        elif atom.symbol != "H":
            raise SystemExit("Error, there are heteroatoms in your structure!")
    carbon_nl = NeighborList(natural_cutoffs(carbon), self_interaction=False, bothways=True)
    carbon_nl.update(carbon)
    con_mat_sparse = carbon_nl.get_connectivity_matrix()
    con_mat = con_mat_sparse.toarray()
    if len(con_mat) %2 != 0:
        raise SystemExit("Error, the molecule has an odd number of carbon atoms!")
    for row in con_mat: 
        if np.sum(row) < 2: raise SystemExit("Error, one carbon is not part of a ring!")
    return con_mat
    
def con_mat_from_SMILES(smiles): # connectivity matrix of carbon structure from SMILES string
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return [] #raise SystemExit("Error, invalid SMILES string!")
    mol = Chem.AddHs(molecule)
    carbon = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
                raise SystemExit("Error, one carbon is not sp2 hybridized!")
            carbon.append(atom.GetIdx())
        elif atom.GetSymbol() != "H":
            raise SystemExit("Error, there are heteroatoms in your structure!")
    con_mat = np.zeros((len(carbon), len(carbon)), dtype=int)
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        if atom1 in carbon and atom2 in carbon:
            con_mat[atom1][atom2] = 1; con_mat[atom2][atom1] = 1
    if len(con_mat) %2 != 0:
        raise SystemExit("Error, the molecule has an odd number of carbon atoms!")
    for row in con_mat: 
        if np.sum(row) < 2: raise SystemExit("Error, one carbon is not part of a ring!")
    return con_mat
    
def con_mat_from_mat(str_mat):
    con_mat = np.array(ast.literal_eval(str_mat))
    if len(con_mat) %2 != 0:
        raise SystemExit("Error, the molecule has an odd number of carbon atoms!")
    for row in con_mat: 
        if np.sum(row) < 2: raise SystemExit("Error, one carbon is not part of a ring!")
        elif np.sum(row) > 3: raise SystemExit("Error, one carbon has too many bonds!")
    return con_mat
