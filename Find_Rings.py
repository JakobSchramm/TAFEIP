def neighbor_dic_from_con_mat(con_mat): # get dictionary of connected atoms from connectivity matrix
    nl_dict = {}
    for i in range(len(con_mat)):
        nl_dict[i] = []
        for j in range(len(con_mat)):
            if con_mat[i, j] == 1:
                nl_dict[i].append(j)
    return nl_dict

def edges_from_nl_dict(nl_dict): # get edges from dictionary
    edges = []
    for i in list(nl_dict):
        for j in nl_dict[i]:
            if sorted([i,j]) not in edges:
                edges.append(sorted([i,j]))
    return edges

def find_all_paths(graph, start, end, path =[]): # find all paths between two points in graph
    path = path + [start]
    if start == end:
        return [path]
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    paths.sort(key=len)
    return paths

def find_possible_circuits(con_mat, needs_to_be_even):
    nl_dict = neighbor_dic_from_con_mat(con_mat) # create neighbor dictionary
    edges = edges_from_nl_dict(nl_dict) # find all existing edges in graph
    
    # find all possible circuits using first step of Hortons algorithm,
    # i.e. for every point (v) and independent edge (wx) combination, 
    # all possible circuits using all possible paths from v to w and v to x are created,
    # and if it is a true circuit and not noticed yet, it is appended to a list and sorted from smallest to largest
    circuits = []
    edges_of_circuit_sorted = []
    for s in list(nl_dict):
        for edge in edges:
            if s not in edge:
                path1 = find_all_paths(nl_dict, s, edge[0])
                path2 = find_all_paths(nl_dict, s, edge[1])
                for p in path1:
                    for q in path2:
                        edges_of_p = []
                        for ind_of_p in range(len(p)-1):
                            edges_of_p.append(sorted([p[ind_of_p], p[ind_of_p+1]]))
                        edges_of_q = []
                        for ind_of_q in range(len(q)-1):
                            edges_of_q.append(sorted([q[ind_of_q], q[ind_of_q+1]]))
                        possible_circuit = p+q[::-1]
                        possible_circuit.pop(-1)
                        possible_edges_of_circuit_sorted = sorted(edges_of_p+edges_of_q+[edge])
                        if len(set(possible_circuit)) == len(possible_circuit):
                            if needs_to_be_even == True:
                                if len(possible_circuit) %2 == 0:
                                    if possible_edges_of_circuit_sorted not in edges_of_circuit_sorted:
                                        circuits.append(possible_circuit)
                                        edges_of_circuit_sorted.append(possible_edges_of_circuit_sorted)
                            else:
                                if possible_edges_of_circuit_sorted not in edges_of_circuit_sorted:
                                    circuits.append(possible_circuit)
                                    edges_of_circuit_sorted.append(possible_edges_of_circuit_sorted)
    circuits.sort(key=len)
    
    # for all circuits, a list of edges for the circuits is created
    edges_of_circuits = []
    for pos1,circuit in enumerate(circuits):
        edges_of_circuits.append([])
        for pos2,v in enumerate(circuit):
            if pos2 < len(circuit)-1:
                edges_of_circuits[pos1].append([circuit[pos2],circuit[pos2+1]])
            else:
                edges_of_circuits[pos1].append([circuit[pos2],circuit[0]])
    
    return edges, circuits, edges_of_circuits
    
    # for all circuits, a list of edges for the circuits is created
    edges_of_circuits = []
    for pos1,circuit in enumerate(circuits):
        edges_of_circuits.append([])
        for pos2,v in enumerate(circuit):
            if pos2 < len(circuit)-1:
                edges_of_circuits[pos1].append([circuit[pos2],circuit[pos2+1]])
            else:
                edges_of_circuits[pos1].append([circuit[pos2],circuit[0]])
    
    return edges, circuits, edges_of_circuits
