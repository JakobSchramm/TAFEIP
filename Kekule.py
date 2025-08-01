import numpy as np

def find_neighbors(atom, mat): # neighbors of atom and types of bonds between atom and neighbor: 1-connected, 10-sibo, 20-dobo
    neighbor = [n for n,x in enumerate(mat[atom]) if x != 0]
    neighbor_type = [x for n,x in enumerate(mat[atom]) if x != 0]
    return neighbor,neighbor_type

def check_valency(atom, mat): # check how many bonds the atom has
    neighbor, neighbor_type = find_neighbors(atom, mat)
    if len(neighbor) == 2: # tertiary atom
        if sum(neighbor_type) == 2: # 2x no bond -> irrelevant
            new_bond_type = 0
        elif sum(neighbor_type) == 11: # 1x sibo & 1x no bond -> dobo
            new_bond_type = 20
        elif sum(neighbor_type) == 21: # 1x dobo & 1x no bond -> sibo
            new_bond_type = 10
        elif sum(neighbor_type) == 30:
            print("Should not happen")
    elif len(neighbor) == 3: # quaternary atom
        if sum(neighbor_type) == 3: # 3x no bond  -> irrelevant
            new_bond_type = 0
        elif sum(neighbor_type) == 12: # 1x sibo & 2x no bond -> choice
            new_bond_type = 100
        elif sum(neighbor_type) == 21: # 2x sibo & 1x no bond -> dobo
            new_bond_type = 20
        elif sum(neighbor_type) == 22: # 1x dobo & 2x no bond -> sibo
            new_bond_type = 10
        elif sum(neighbor_type) == 31: # 1x sibo & 1x dobo & 1x no bond -> sibo
            new_bond_type = 10
        elif sum(neighbor_type) == 40:
            print("Should not happen")
    return new_bond_type

def check_choice(atom1, atom2, mat): # check between two atoms, what bonds they have and what they thus want
    new_bond_type1 = check_valency(atom1, mat)
    new_bond_type2 = check_valency(atom2, mat)
    if new_bond_type1 + new_bond_type2 == 0 or new_bond_type1 + new_bond_type2 == 100 or new_bond_type1 + new_bond_type2 == 200: # both are irrelevant, one has a choice and the other is irrelevent, or both have a choice
        new_bond_type = "choice"
    elif new_bond_type1 + new_bond_type2 == 30: # one atom needs a different bond then the other -> discard this structure later
        new_bond_type = "Error"
    elif new_bond_type1 + new_bond_type2 == 110 or new_bond_type1 + new_bond_type2 == 120: # one has choice and the other needs a bond -> give the bond
        new_bond_type = new_bond_type1 + new_bond_type2 - 100
    elif new_bond_type1 == 0 or new_bond_type2 == 0: # one is irrelevant and the other needs a bond -> give the bond
        new_bond_type = new_bond_type1 + new_bond_type2
    else: # remaining case (10==10 or 20==20): both need the same bond -> give this bond
        new_bond_type = new_bond_type1
    return new_bond_type

def add_sibo(atom1, atom2, mat): # add single bond if they are connected
    if mat[atom1][atom2] != 0: mat[atom1][atom2] = 10; mat[atom2][atom1] = 10
    else: raise SystemExit("Error, atoms that should be connected are not connected!") # sanity check
    return mat

def add_dobo(atom1, atom2, mat): # add double bond if they are connected
    if mat[atom1][atom2] != 0: mat[atom1][atom2] = 20; mat[atom2][atom1] = 20
    else: raise SystemExit("Error, atoms that should be connected are not connected!") # sanity check
    return mat

def add_one_bond(atom1, Kekule_mat): # add one bond to current atom or if saturated go to next
    neighbor, neighbor_type = find_neighbors(atom1, Kekule_mat) # find all neighbors
    if sum(neighbor_type) != 30 and sum(neighbor_type) != 40: # if the current atom is not saturated
        atom2 = neighbor[neighbor_type.index(1)] # find neighbor to which no bond is yet
        new_bond_type = check_choice(atom1, atom2, Kekule_mat) # check which bond should be there
        if new_bond_type == "Error": # complication in assignment: one atom needs a different bond then the other -> discard this structure
            list_Kekule_mat = []
        elif new_bond_type == 20: # add dobo
            Kekule_mat = add_dobo(atom1, atom2, Kekule_mat)
            list_Kekule_mat = [Kekule_mat]
        elif new_bond_type == 10: # add sibo
            Kekule_mat = add_sibo(atom1, atom2, Kekule_mat)
            list_Kekule_mat = [Kekule_mat]
        else: # if there is a choice, create 2 copies and a single and add a double bond, respectively
            Kekule_mat10 = np.copy(Kekule_mat)
            Kekule_mat20 = np.copy(Kekule_mat)
            Kekule_mat10 = add_sibo(atom1, atom2, Kekule_mat10)
            Kekule_mat20 = add_dobo(atom1, atom2, Kekule_mat20)
            list_Kekule_mat = [Kekule_mat10, Kekule_mat20]
    else: # if current atom is saturated
        for neighbor_index, neighbor_atom in enumerate(neighbor): # search for new neighbor that is still not saturated
            new_neighbor, new_neighbor_type = find_neighbors(neighbor_atom, Kekule_mat)
            list_Kekule_mat = [Kekule_mat]
            if 1 in new_neighbor_type: atom1 = neighbor_atom; break # choose new atom and stop if one finds one
            else: atom1 = "end" # no new atoms can be found
    return atom1, list_Kekule_mat

def manipulate_Kekule(atom1, list_Kekule_mat): # manipulate all structures in one step
    new_atom1 = "new is needed"
    for Kekule_num, Kekule_mat in enumerate(list_Kekule_mat.copy()): # for all current Kekule structures in list
        if set(map(tuple,Kekule_mat)) != set(map(tuple,[["failed"]])): # only if the current matrix is not one that had a complication in assignment
            new_atom1, new_Kekule_mat =  add_one_bond(atom1, Kekule_mat) # manipulate accordingly and get new atom
            if len(new_Kekule_mat) == 2: # if there was a choice, two matrices are given out.
                list_Kekule_mat[Kekule_num] = new_Kekule_mat[0] # replace the old one if the first choice
                list_Kekule_mat.append(new_Kekule_mat[1]) # append the second choice
            elif len(new_Kekule_mat) == 1: # if there was no choice
                list_Kekule_mat[Kekule_num] = new_Kekule_mat[0] # replace the old one with new one
            elif len(new_Kekule_mat) == 0: # if there was a complication in assignment
                list_Kekule_mat[Kekule_num] = [["failed"]] # add failed comment
    if new_atom1 != "new is needed": # check if a new atom is defined where one can continue
    	atom1 = new_atom1
    else:
        atom1 = "end"
    return atom1, list_Kekule_mat

def find_Kekule(atom1, con_mat): # find Kekule in structure - only in substructure can happen!
    list_Kekule_mat = [np.copy(con_mat)]
    while atom1 != "end": # as long as atom is not saturated
        atom1, list_Kekule_mat = manipulate_Kekule(atom1, list_Kekule_mat) # manipulate Kekule structure accordingly
    final_list_Kekule_mat = [mat for mat in list_Kekule_mat if set(map(tuple,mat)) != set(map(tuple,[["failed"]]))] # remove all failed attempts
    return final_list_Kekule_mat

def check_if_done(list_Kekule_mat): # check if there are atoms with no bond assigned (1) left
    sublist_Kekule_mat = []
    check_again = False
    for Kekule_mat in list_Kekule_mat.copy():
        missing_list = []
        for row_ind,row in enumerate(Kekule_mat):
            for element_ind,element in enumerate(row):
                if element == 1 and element_ind not in missing_list: # if there is only connection
                    missing_list.append(element_ind) # add to the list of missing atoms
        if missing_list != []: # if this list is not empty
            for atom in missing_list:
                neighbor, neighbor_type = find_neighbors(atom, Kekule_mat)
                if len(neighbor) == 2: new_atom = atom; break # find first atom that is tertiary, stop searching and continue with that atom
                else: new_atom = atom # if no tertiary atom is possible, continue with last quarternary atom
            sublist_Kekule_mat.extend(find_Kekule(new_atom, Kekule_mat)) # find substructures
            check_again = True # check later again if there are atoms with no bond assigned
        else: # if missinglis is empty - there are no unassigned bonds left
            sublist_Kekule_mat = list_Kekule_mat # safe the old list
            check_again = False # no need to check again
    return sublist_Kekule_mat, check_again

def find_all_Kekule(con_mat): # find all Kekule structures given a connectivity matrix
    for atom_ind, row in enumerate(con_mat): # look for tertiary starting atom
        if sum(row) == 2: start_atom = atom_ind; break
    list_Kekule_mat = find_Kekule(start_atom, con_mat) # find Kekule structure starting from starting atom
    check_again = True
    while check_again: # check if there are atoms with no bond assigned left
        list_Kekule_mat, check_again = check_if_done(list_Kekule_mat)
    if list_Kekule_mat == []:
        raise SystemExit("Error, the molecule has no Kekule structure!")
    return list_Kekule_mat
