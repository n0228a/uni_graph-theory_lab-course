# A DRF returns the symmetric difference between educt and product graphs == Reaction graph
from hashlib import blake2b
from synkit.IO import rsmi_to_graph
from phi_transformation import phi_edge_graph, phi_shortest_path_graph, phi_vertex_graph


def vertex_drf(rsmi:str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)
  
  # Calculate vertex sets for educt and product graph
  educt_vertices = phi_vertex_graph(educt_graph)
  product_vertices = phi_vertex_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_vertices).symmetric_difference(set(product_vertices))

  return educt_vertices, product_vertices, sym_diff

def edge_drf(rsmi:str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)
  
  # Calculate edge sets for educt and product graph
  educt_edges = phi_edge_graph(educt_graph)
  product_edges = phi_edge_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_edges).symmetric_difference(set(product_edges))

  return educt_edges, product_edges, sym_diff

def shortest_path_drf(rsmi:str):
  educt_graph, product_graph = rsmi_to_graph(rsmi)
  
  # Calculate shortest path sets for educt and product graph
  educt_paths = phi_shortest_path_graph(educt_graph)
  product_paths = phi_shortest_path_graph(product_graph)

  # Calculate symmetric difference
  sym_diff = set(educt_paths).symmetric_difference(set(product_paths))

  return educt_paths, product_paths, sym_diff

def summary(before, after, smy_diff):
  print(f"Before: {before}")
  print(f"After: {after}")
  print(f"Symmetric Difference: {smy_diff}")