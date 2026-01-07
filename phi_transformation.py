from synkit.IO import rsmi_to_graph
from _blake2 import blake2b
from networkx.algorithms import all_pairs_shortest_path

hasher = blake2b(digest_size=32)

# Transforms the vertices of a SMILES reaction's educt graph into their hashed representation
# Returns a set of hashed, concated vertex labels for the educt graph of the passed SMILES reaction
# E.g. {hashed_label1, hashed_label2hashed_label3, ...}
def phi_vertex(rsmi: str = "[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1"):
    graph, _ = rsmi_to_graph(rsmi)
    vertex_set = {}
    vertex_labels = set()
    for n, d in graph.nodes(data=True):
        if n not in vertex_set:
            vertex_set[n] = d['element']
        else:
            # Add d['element'] to existing entry
            vertex_set[n] += d['element']
    
    # Hash the labels
    for value in vertex_set.values():
        hasher.update(value.encode('utf-8'))
        hashed = hasher.digest().hex()  # Get the actual hash
        vertex_labels.append(hashed)  # Append the hashed value
  
    return vertex_labels

# Transforms the vertices of a SMILES reaction's educt graph into their hashed representation
# Returns a dict of hashed vertices with their hashed labels for the educt graph of the passed SMILES reaction
# E.g. {hashed_node1: {hashed_label1, hashed_label2}, hashed_node2: {hashed_label3}, ...}
def phi_vertex_dict(rsmi: str = "[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1"):
    graph, _ = rsmi_to_graph(rsmi)
    vertex_set = dict()
    vertex_labels = dict()
    for n, d in graph.nodes(data=True):
        if n not in vertex_set:
            vertex_set[n] = set()
            vertex_set[n].add(d['element'])
        else:
            vertex_set[n].add(d['element'])
    
    for key, value_set in vertex_set.items():
        # Hash the key
        hasher.update(str(key).encode('utf-8'))
        hashed_key = hasher.digest().hex()  # Get the actual hash
        # transform all items in value_set to their hashes
        hashed_values = set()
        for val in value_set:
            hasher.update(val.encode('utf-8'))
            hashed_val = hasher.digest().hex()  # Get the actual hash
            hashed_values.add(hashed_val)

        # Sort hashed_values to ensure consistent representation
        hashed_values = frozenset(sorted(hashed_values))
        vertex_labels[hashed_key] = hashed_values  # Append the hashed value
  
    return vertex_labels

# Transform the edges of a SMILES reaction's educt graph into their hashed representation
# The plain representation of an edge is "uorderv" where u and v are the node indices and order is the bond order
# Returns a set of hashed edges for the educt graph of the passed SMILES reaction
# E.g. {hashed_edge1, hashed_edge2, ...} wheras hashed_edge is, for example, 131.02 with "13" and "2" being the node indices and "1.0" the bond order
def phi_edge(rsmi: str = "[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1"):
    graph, _ = rsmi_to_graph(rsmi)
    edge_set = list()

    for u, v, d in graph.edges(data=True):
        # Sort u and v to ensure consistent representation
        u, v = sorted([u, v])
        result = f"{u}{d['order']}{v}"
        # Hash edge representation
        hasher.update(result.encode('utf-8'))
        hashed = hasher.digest().hex()  # Get the actual hash
        edge_set.append(hashed)  # Append the hashed value
  
    return edge_set

# Computes all shortest paths for a SMILES reaction's educt graph and returns their hashed representation
# Returns a set of hashed shortest paths for the educt graph of the passed SMILES reaction
# Default for testing purposes
def phi_shortest_path(rmsi: str = "[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1"):
    graph, _ = rsmi_to_graph(rmsi)
    graph = graph.to_undirected() # Should naturally be undirected, but just to be sure
    node_to_label = {n: d["element"] for n, d in graph.nodes(data=True)}
    paths = dict(all_pairs_shortest_path(graph))
    paths_set = set()
    paths_set_hased = set()

    # For each path, convert to string representation
    for source, target_dict in paths.items():
        for target, path in target_dict.items():
            # Convert the path to label concatenation using node_to_label
            label_path = ''.join(node_to_label[n] for n in path)
            # Add source and target to label_path
            label_path = f"{node_to_label[source]}{label_path}{node_to_label[target]}"
            paths_set.add(label_path)
            # Hash label_path
            hasher.update(label_path.encode('utf-8'))
            hashed = hasher.digest().hex()  # Get the actual hash
            paths_set_hased.add(hashed)  # Append the hashed value

    return paths_set_hased  

############ SAME METHODS WITH GRAPH AS ARG ############

# Transforms the vertices of a passed graph into their hashed representation
# Returns a set of hashed, concated vertex labels for the educt graph of the passed SMILES reaction
# E.g. {hashed_label1, hashed_label2hashed_label3, ...}
def phi_vertex_graph(graph):
    vertex_set = {}
    vertex_labels = list()
    for n, d in graph.nodes(data=True):
        if n not in vertex_set:
            vertex_set[n] = d['element']
        else:
            # Add d['element'] to existing entry
            vertex_set[n] += d['element']
    
    # Hash the labels
    for value in vertex_set.values():
        hasher.update(value.encode('utf-8'))
        hashed = hasher.digest().hex()  # Get the actual hash
        vertex_labels.append(hashed)  # Append the hashed value
  
    return vertex_labels

# Transforms the vertices of a passed graph into their hashed representation
# Returns a dict of hashed vertices with their hashed labels as set
# E.g. {hashed_node1: {hashed_label1, hashed_label2}, hashed_node2: {hashed_label3}, ...}
def phi_vertex_dict_graph(graph):
    vertex_set = dict()
    vertex_labels = dict()
    for n, d in graph.nodes(data=True):
        if n not in vertex_set:
            vertex_set[n] = set()
            vertex_set[n].add(d['element'])
        else:
            vertex_set[n].add(d['element'])
    
    for key, value_set in vertex_set.items():
        # Hash the key
        hasher.update(str(key).encode('utf-8'))
        hashed_key = hasher.digest().hex()  # Get the actual hash
        # transform all items in value_set to their hashes
        hashed_values = set()
        for val in value_set:
            hasher.update(val.encode('utf-8'))
            hashed_val = hasher.digest().hex()  # Get the actual hash
            hashed_values.add(hashed_val)

        # Sort hashed_values to ensure consistent representation
        hashed_values = frozenset(sorted(hashed_values))
        vertex_labels[hashed_key] = hashed_values  # Append the hashed value
  
    return vertex_labels

# Transform the edges of a passed graph into their hashed representation
# The plain representation of an edge is "uorderv" where u and v are the node indices and order is the bond order
# Returns a set of hashed edges for the educt graph of the passed SMILES reaction
# E.g. {hashed_edge1, hashed_edge2, ...} wheras hashed_edge is, for example, 131.02 with "13" and "2" being the node indices and "1.0" the bond order
def phi_edge_graph(graph):
    edge_set = set()

    for u, v, d in graph.edges(data=True):
        # Sort u and v to ensure consistent representation
        u, v = sorted([u, v])
        result = f"{u}{d['order']}{v}"
        # Hash edge representation
        hasher.update(result.encode('utf-8'))
        hashed = hasher.digest().hex()  # Get the actual hash
        edge_set.add(hashed)  # Add the hashed value
  
    return edge_set

# Computes all shortest paths for the passed graph and returns their hashed representation
# Returns a set of hashed shortest paths for the passed graph
# Default for testing purposes
def phi_shortest_path_graph(graph):
    graph = graph.to_undirected() # Should naturally be undirected, but just to be sure
    node_to_label = {n: d["element"] for n, d in graph.nodes(data=True)}
    paths = dict(all_pairs_shortest_path(graph))
    paths_set = set()
    paths_set_hased = set()
    # For each path, convert to string representation
    for source, target_dict in paths.items():
        for target, path in target_dict.items():
            # Convert the vertice-path to the label-path
            label_path = ''.join(node_to_label[n] for n in path)
            # Extend label_path by adding source and target
            label_path = f"{node_to_label[source]}{label_path}{node_to_label[target]}"
            hasher.update(label_path.encode('utf-8'))
            hashed = hasher.digest().hex()  # Get the actual hash
            paths_set_hased.add(hashed)  # Append the hashed value

    print(paths_set)
  
    return paths_set_hased