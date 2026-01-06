from synkit.IO import rsmi_to_graph
from _blake2 import blake2b
from networkx.algorithms import all_pairs_shortest_path

hasher = blake2b(digest_size=16)

def phi_vertex(rsmi: str):
  graph, _ = rsmi_to_graph(rsmi)
  vertex_set = {}
  vertex_labels = list()
  for n, d in graph.nodes(data=True):
    if n not in vertex_set:
      vertex_set[n] = d['element']
    else:
      # Add d['element'] to existing entry
      vertex_set[n] += d['element']
  # Return hashed representation of vertex_set
  
  for key, value in vertex_set.items():
    # Hash the labels
    hasher.update(value.encode('utf-8'))
    vertex_labels.append(value)

  return vertex_labels


def k_vertex(rsmi: str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)

  count = 0

  for n, d in educt_graph.nodes(data=True):
    for np, dp in product_graph.nodes(data=True):
      if d['element'] == dp['element']:
        count += 1

  return count

def phi_edge(rsmi: str):
  graph, _ = rsmi_to_graph(rsmi)

  edge_set = list()

  for u, v, d in graph.edges(data=True):
    result = f"{u}{d['order']}{v}"
    # Hash edge representation
    hasher.update(result.encode('utf-8'))
    edge_set.append(result)
  
  return edge_set

def phi_shortest_path(rmsi:str):
  graph, _ = rsmi_to_graph(rmsi)
  # One path is the concatenation of all labels along the path from node A to node B

  for n, d in graph.nodes(data=True):
    node_to_label = {n: d["element"] for n, d in graph.nodes(data=True)}

  paths = dict(all_pairs_shortest_path(graph))
  paths_set = list()

  # For each path, convert to string representation
  for source, target_dict in paths.items():
    for target, path in target_dict.items():
      # Convert the path to label concatenation using node_to_label
      label_path = ''.join(node_to_label[n] for n in path)
      # Add source and target to label_path
      label_path = f"{node_to_label[source]}{label_path}{node_to_label[target]}"
      # Hash label_path
      hasher.update(label_path.encode('utf-8'))
      paths_set.append(label_path)
  
  return paths_set

# [CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1