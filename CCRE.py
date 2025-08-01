import numpy as np
from Kekule import find_all_Kekule
from Find_Rings import find_possible_circuits
import itertools
    
def check_alternating(bonds): # check if bonds are alternating, i.e. if two neighboring bonds do not have the same type
    for i in range(0, len(bonds)-1):
        if bonds[i] == bonds[i+1]:
            return False
    if bonds[i] == bonds[-1]:
        return False
    return True
    
def find_conj_circuits(Kekule_mat, edges_of_circuits): # check for every circuit if it is alternating
    conjugated_circuits = []
    for circuit in edges_of_circuits:
        bonds = []
        for atom1, atom2 in circuit: # append bond type of all bonds in circuit to list
            bonds.append(Kekule_mat[atom1][atom2])
        alternating = check_alternating(bonds) # check
        if alternating: conjugated_circuits.append(circuit) # if alternating, add to conjugated circuits
    return conjugated_circuits
    
def check_linear_dependence(edges, conjugated_circuits): # check if ANY of the longest conjugated circuits is a linear combination of all other conjugated circuits
    conjugated_circuits.sort(key=len, reverse=True) # sort in reverse after the length
    
    # sort the edges in conjugated circuits
    conjugated_circuits_sorted = []
    for c in conjugated_circuits:
        c_sorted = []
        for e in c:
            c_sorted.append(sorted(e))
        conjugated_circuits_sorted.append(sorted(c_sorted))
    
    # create a vector for all edges that contains 1 if edge in conjugated circuit, else 0 (for all conjugated circuits)
    all_edge_vectors = []
    for conjucated_circuit in conjugated_circuits_sorted:
        edge_vector = [0]*len(edges)
        for e_ind,edge in enumerate(edges):
            if edge in conjucated_circuit:
                edge_vector[e_ind]=1
        all_edge_vectors.append(np.array(edge_vector))
    
    # check if ANY of the long conjugated circuits is a linear combination of smaller conjugated circuits
    linear_combs = list(map(list, itertools.product([1, 0, -1], repeat=len(all_edge_vectors)-1))) # create all combination factors for other circuits
    remove_linear_comb_list = []
    for length_index,edge_vector in enumerate(all_edge_vectors): # iterate over all circuits
        edge_vector_to_test = edge_vector # get vector to test
        remaining_edge_vectors = all_edge_vectors[length_index+1:] # get remaining smaller vectors
        for comb in linear_combs: # iterate over all possible combination factors
            sum = np.array([0]*len(edges))
            for edge_ind, edge_vector in enumerate(remaining_edge_vectors):
                sum += edge_vector * comb[edge_ind] # create combination iteratively
            if list(sum) == list(edge_vector_to_test): # check if combination is equal to vector to test
                remove_linear_comb_list.append(length_index); break # if yes, remove the circuit. On found combination is enough evidence
    new_conjugated_circuits = sorted([conj_circuit for circuit_ind, conj_circuit in enumerate(conjugated_circuits) if circuit_ind not in remove_linear_comb_list], key=len) # if it is linear combination, remove conjugated circuits
    return new_conjugated_circuits
    
def calc_resonance_vector(list_Kekule, edges_of_circuits, edges, resonance_rings):
    total_resonance_vector = np.array([0]*len(resonance_rings))
    max_number_of_conjugated_circuits = 0
    number_of_circuits = []
    for Kekule_mat in list_Kekule:
        conjugated_circuits = find_conj_circuits(Kekule_mat, edges_of_circuits) # find all conjugated circuits
        if len(conjugated_circuits) > max_number_of_conjugated_circuits: # save the maximum number of conjugated circuits including linear combinations for later sanity check 1
            max_number_of_conjugated_circuits = len(conjugated_circuits)
        conjugated_circuits = check_linear_dependence(edges, conjugated_circuits)
        number_of_circuits.append(len(conjugated_circuits)) # save the number of linear independent circuits of each Kekule structure for later sanity check 2
        resonance_list = [len(c) for c in conjugated_circuits] # get the lengths of the conjugated circuits
        resonance_vector = np.array([0]*len(resonance_rings)) # create vector
        for r_ind,r in enumerate(resonance_rings):
            resonance_vector[r_ind] = resonance_list.count(r) # add number of circuits per length
        total_resonance_vector += resonance_vector # sum all up
    if max_number_of_conjugated_circuits != len(list_Kekule)-1: # sanity check 1
        print("### WARNING ###")
        print("The maximum number of conjugated circuits is not 1 less then the number of Kekule structures.")
        print("However, this might be an artefect as there might be linear combinations of smaller circuits missing.")
    if len(set(number_of_circuits)) != 1: # sanity check 2
        raise SystemExit("Error, the number of linear independent circuits is different between different Kekule structures.This should never happen! (The script might have failed to exclude a linear combination.)")
    return total_resonance_vector
    
def R_or_Q(PE): # Calculation of resonance energy of the ciruit using the HMO parametrization
    sum = 0
    for i in range(1,int(PE/4)+1):
        sum += np.cos(2 * np.pi * i / PE)
    E_del = 4 + 8 * sum - PE
    E_inf_per_PE = 4 / np.pi - 1
    R_or_Q = E_del - E_inf_per_PE * PE
    return R_or_Q 

def calc_RE(con_mat): # Calculate lists of all possible circuit sizes (up to number of PE as this is largest possible cycle) and resonance energies
    resonance_rings = []
    resonance_energy = []
    for i in range(4, len(con_mat)+1, 2):
        resonance_rings.append(i)
        resonance_energy.append(R_or_Q(i))
    return resonance_rings, resonance_energy
    
def calc_CCRE(con_mat): # calculate CCRE from vector of circuit lengths and print CC expression
    resonance_rings, resonance_energy = calc_RE(con_mat)
    list_Kekule = find_all_Kekule(con_mat)
    edges, circuits, edges_of_circuits = find_possible_circuits(con_mat, True)
    if circuits == []: # sanity check
        raise SystemExit("Error, no rings were detected!")
    total_resonance_vector = calc_resonance_vector(list_Kekule, edges_of_circuits, edges, resonance_rings)
    
    letters = [f"R_{int((ind+1)/2)}" if ind%2 !=0 else f"Q_{int((ind+2)/2)}" for ind,val in enumerate(total_resonance_vector)]
    print("Expression for CCRE: ( "+" + ".join([f"{num} {letters[i]}" for i, num in enumerate(total_resonance_vector) if num > 0])+f" ) / {len(list_Kekule)}")
    
    CCRE = np.sum(total_resonance_vector * np.array(resonance_energy)) / len(list_Kekule)
    CCREPE = CCRE / len(con_mat)
    return CCRE, CCREPE #, total_resonance_vector, len(list_Kekule)
