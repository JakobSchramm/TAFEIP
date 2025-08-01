import numpy as np
import sympy as sp
from sympy import Matrix, simplify, nroots, Poly
from sympy.abc import x
from Find_Rings import find_possible_circuits

def find_quaternary(con_mat): # find all quarternary carbon atoms and print to list
    quart_list = []
    for atom,row in enumerate(con_mat):
        if np.sum(row) == 3:
            quart_list.append(atom)
    return quart_list
    
def find_ring_and_or_quart(con_mat): # find all carbon atoms that are part of a ring, quarternary, or both
    edges, circuits, edges_of_circuits = find_possible_circuits(con_mat, False) # find rings
    quart_list = find_quaternary(con_mat) # find quarternary carbons
    ring_and_quart = []
    only_ring = []
    only_quart = []
    for atom in range(len(con_mat)):
        if any(atom in circuit for circuit in circuits) and atom in quart_list: # if both
            ring_and_quart.append(atom)
        elif any(atom in circuit for circuit in circuits) and atom not in quart_list: # if only ring
            only_ring.append(atom)
        elif atom in quart_list and any(atom not in circuit for circuit in circuits): # if only quarternary
            only_quart.append(atom)
    return ring_and_quart, only_ring, only_quart, edges
   
def remove_edge(atom1, atom2, mat): # remove the edge between atom1 and atom2
    mat_copy = np.copy(mat)
    if mat_copy[atom1][atom2] == 1 and mat_copy[atom2][atom1] == 1: # sanity check
        mat_copy[atom1][atom2] = 0; mat_copy[atom2][atom1] = 0
    else:
        raise SystemExit("Error, the bond that should be removed is not a bond!")
    return mat_copy

def remove_edge_with_vert(atom1, atom2, mat): # remove the vertices of atom1 and atom2 and all associated edges. In fact, we are not removing the vertices, so an additional x appears later in the characteristic polynomial!
    mat_copy = np.copy(mat)
    if mat_copy[atom1][atom2] == 1 and mat_copy[atom2][atom1] == 1: # sanity check
        neighbors_atom1 = [n for n,x in enumerate(mat_copy[atom1]) if x != 0] # find neighbors of atom1
        neighbors_atom2 = [n for n,x in enumerate(mat_copy[atom2]) if x != 0] # find neighbors of atom2
        for neighbor_atom in neighbors_atom1:
            mat_copy[neighbor_atom][atom1] = 0; mat_copy[atom1][neighbor_atom] = 0 # remove edges of atom1 to neighbors
        for neighbor_atom in neighbors_atom2:
            mat_copy[neighbor_atom][atom2] = 0; mat_copy[atom2][neighbor_atom] = 0 # remove edges of atom2 to neighbors
        mat_copy[atom1][atom2] = 0; mat_copy[atom2][atom1] = 0 # remove edge between atom1 and atom2
    else:
        raise SystemExit("Error, the bond that should be removed is not a bond!")
    return mat_copy
    
def find_acyclic_mats(con_mat): # find all matrices of acyclic fragments
    acyclic_mats = [[1,np.copy(con_mat),[]]] # initialize with connectivity matrix: [sign, matrix, removed_atoms]
    check_again = [True]
    while set(check_again) != set([False]): # do as long as we are manipulating the list of all matrices of acyclic fragments
        check_again = []
        new_acyclic_mats = []
        for mat_num,(sign,mat,removed_atoms) in enumerate(acyclic_mats): # loop over all matrices of acyclic fragments
            ring_and_quart, only_ring, only_quart, edges = find_ring_and_or_quart(mat) # find the interesting atoms
            if ring_and_quart != []: # if there is an atom that is part of a ring and quarternary, start there
                atom1 = ring_and_quart[0]
                manipulate_mat = True # we need to manipulate the original matrix as it is not linear yet 
            elif only_ring != []: # else, if there is an atom that is only part of a ring, start there
                atom1 = only_ring[0]
                manipulate_mat = True
            elif only_quart != []: # else, if there is an atom that is only quarternary, start there
                atom1 = only_quart[0]
                manipulate_mat = True
            else:  # else we don't need to manipulate the original matrix
                manipulate_mat = False
            if manipulate_mat: # if we need to manipulate the matrix
                for potential_neighbor in edges: # find a neighbor in all edges
                    if atom1 == potential_neighbor[0]:
                        atom2 = potential_neighbor[1]
                        break # stop if first neighbor is found
                    elif atom1 == potential_neighbor[1]:
                        atom2 = potential_neighbor[0]
                        break # same
                new_acyclic_mats.append([sign, remove_edge(atom1, atom2, mat), removed_atoms]) # remove the edge, keep the sign and previously removed atoms
                new_acyclic_mats.append([-sign, remove_edge_with_vert(atom1, atom2, mat), removed_atoms+[atom1,atom2]]) # remove the vertices (atoms) and all associated edges, change the sign and add removed atoms to previously removed atoms
                check_again.append(True) # we need to do this again as new matrices are appended
            else:
                new_acyclic_mats.append([sign,mat,removed_atoms]) # if we did not manipulate the matrix, we just append it as it is
                check_again.append(False) # we don't need to check this matrix again
        acyclic_mats = new_acyclic_mats # the matrices of acyclic fragments can be initialized with the new matrices of acyclic fragments
    return acyclic_mats
    
def find_roots(char_poly):
    roots = []
    for root in nroots(char_poly, n=8, maxsteps=200):
        if type(root) == sp.core.add.Add:
            if root.args[1].args[0] < 1e-3:
                print("### WARNING ###")
                print(f"A root of a polynomial had an imaginary part: {root}")
                print("However, as it is small it is neglected.")
                roots.append(root.args[0])
            else:
                raise SystemExit(f"Error, there is a root that as a large imaginary part: {root}")
        else:
            roots.append(root)
    roots.sort(reverse = True)
    return roots

def calc_TRE(con_mat):
    acyclic_mats = find_acyclic_mats(con_mat)
    
    mol_char_poly = Matrix(con_mat).charpoly(x) # calculate the characteristic polynomial of the molecule
    print("Characteristic Polynomial of Molecule:         ", mol_char_poly.as_expr())
    roots_mol = find_roots(mol_char_poly) # find its roots
    
    acyclic_char_poly = 0
    for sign,mat,removed_atoms in acyclic_mats: # calculating the characteristic polynomial of the acyclic reference by adding the polynomials of all linear fragments. Careful, 1 removed atom results in an factor of x, so we need to divide by x^(number of removed atoms)
        acyclic_char_poly += sign * Matrix(mat).charpoly(x) / (x**(len(removed_atoms)))
    print("Characteristic Polynomial of Acyclic Reference:", simplify(acyclic_char_poly).as_expr())
    roots_acyclic = find_roots(acyclic_char_poly) # find its roots
    
    if len(roots_mol) != len(roots_acyclic): # sanity check
        raise SystemExit("Error, the polynomials have a different degree!")
    TRE = 0
    for j in range(int(len(roots_mol)/2)): # calc TRE 
        TRE += 2 * (roots_mol[j] - roots_acyclic[j])
    TREPE = TRE/len(con_mat) # and TREPE
    return TRE, TREPE #, mol_char_poly.all_coeffs(), Poly(simplify(acyclic_char_poly)).all_coeffs()
