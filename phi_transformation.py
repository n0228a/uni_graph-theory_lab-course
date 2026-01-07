from synkit.IO import rsmi_to_graph
from _blake2 import blake2b
from networkx.algorithms import all_pairs_shortest_path
#from multiset import *

hasher = blake2b(digest_size=32)

def phi_vertex(rsmi: str = "[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1"):
    graph, _ = rsmi_to_graph(rsmi)
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
        for val in value_set:
            hasher.update(val.encode('utf-8'))
            hashed_val = hasher.digest().hex()  # Get the actual hash
            value_set.remove(val)
            value_set.add(hashed_val)

        vertex_labels[hashed_key] = value_set  # Append the hashed value
  
    return vertex_labels

def k_vertex(rsmi: str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)

  count = 0

  for n, d in educt_graph.nodes(data=True):
    for np, dp in product_graph.nodes(data=True):
      if d['element'] == dp['element']:
        count += 1

  return count

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
  
    print(paths_set)

    return paths_set_hased  

# [CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1


### SAME METHODS WITH GRAPH AS ARG

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

# Trasforms the edges of a graph into their hashed representation
# Returns a set of hashed edges of the passed graph
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

# Computes all shortest paths within a graph and returns their hashed representation
# Returns a set of hashed shortest paths of the passed graph
def phi_shortest_path_graph(graph):
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
            # Hash label_path
            hasher.update(label_path.encode('utf-8'))
            hashed = hasher.digest().hex()  # Get the actual hash
            paths_set_hased.add(hashed)  # Append the hashed value

    print(paths_set)
  
    return paths_set_hased
# [CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1