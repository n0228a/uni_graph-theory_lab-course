# A DRF returns the symmetric difference between educt and product graphs == Reaction graph
from hashlib import blake2b
from synkit.IO import rsmi_to_graph
from phi_transformation import phi_edge_graph, phi_shortest_path_graph, phi_vertex_dict_graph

def calculate_symmetric_difference_off_dict(dict1, dict2):
  # Calculate for each key in both dicts the symmetric difference of their value sets
  symm_diff = dict()
  all_keys = set(dict1.keys()).union(set(dict2.keys()))
  for key in all_keys:
    set1 = set(dict1.get(key, []))
    set2 = set(dict2.get(key, []))
    symm_diff[key] = set1.symmetric_difference(set2)
  return symm_diff

def vertex_drf(rsmi:str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)
  
  # Calculate vertex sets for educt and product graph
  educt_vertices = phi_vertex_dict_graph(educt_graph)
  product_vertices = phi_vertex_dict_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = calculate_symmetric_difference_off_dict(educt_vertices, product_vertices)

  return sorted(educt_vertices), sorted(product_vertices), sorted(sym_diff)

def edge_drf(rsmi:str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)
  
  # Calculate edge sets for educt and product graph
  educt_edges = phi_edge_graph(educt_graph)
  product_edges = phi_edge_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_edges).symmetric_difference(set(product_edges))

  return sorted(educt_edges), sorted(product_edges), sorted(sym_diff)

def shortest_path_drf(rsmi:str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)
  
  # Calculate shortest path sets for educt and product graph
  educt_paths = phi_shortest_path_graph(educt_graph)
  product_paths = phi_shortest_path_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_paths).symmetric_difference(set(product_paths))

  return sorted(educt_paths), sorted(product_paths), sorted(sym_diff)

######### Same methods with graphs passed directly #########

def vertex_drf_graph(educt_graph, product_graph):
  # Calculate vertex sets for educt and product graph
  educt_vertices = phi_vertex_dict_graph(educt_graph)
  product_vertices = phi_vertex_dict_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = calculate_symmetric_difference_off_dict(educt_vertices, product_vertices)

  return sorted(educt_vertices), sorted(product_vertices), sorted(sym_diff)

def edge_drf_graph(educt_graph, product_graph):
  # Calculate edge sets for educt and product graph
  educt_edges = phi_edge_graph(educt_graph)
  product_edges = phi_edge_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_edges).symmetric_difference(set(product_edges))

  return sorted(educt_edges), sorted(product_edges), sorted(sym_diff)

def shortest_path_drf_graph(educt_graph, product_graph):
  # Calculate shortest path sets for educt and product graph
  educt_paths = phi_shortest_path_graph(educt_graph)
  product_paths = phi_shortest_path_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_paths).symmetric_difference(set(product_paths))

  return sorted(educt_paths), sorted(product_paths), sorted(sym_diff)

