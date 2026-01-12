from phi_transformation import *
from networkx import Graph
from synkit.IO import rsmi_to_graph
from multiset import *

def wfl(graph, h_max, setting):
  # Initial features
  ## Call wfl_vertex, wfl_edge or wfl_shortest_path based on setting

  None

def wfl_vertex(graph, h_max):
  feature_set:dict = phi_vertex_dict_graph_no_hash(graph)

  for h in range(h_max):
    for node in feature_set.keys():
      neighbors = get_nth_neighbors(graph, node, h) # Calculates the neighbors at h distance, return the nodes (int)
      neighbor_labels = list()
      for neighbor in neighbors:
        neighbor_labels.extend(feature_set[neighbor])

      # Append labels to current nodes' feature
      current_labels_of_node = list(feature_set[node])
      current_labels_of_node.extend(neighbor_labels)
      # Sort the labels to ensure consistent ordering
      current_labels_of_node = sorted(current_labels_of_node)

      # Update the feature set
      feature_set[node] = current_labels_of_node
      
  print ("Final Feature Set:", feature_set)

  # TODO: Problem: We've to do wfl for each label separately. Now, we do it simultaneously which is wrong.

  return feature_set

def wfl_vertex_for_node(feature_set, h_max, node):


  # Iterate h_max times
  ## For each node in feature_set, get the set and check for each for neighbors, don't forget to sort
  ## Extend the feature set by appending neighbor features
  None

def get_test_graph():
  g, _ = rsmi_to_graph("[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1")
  return g

def get_nth_neighbors(graph, node, n):
  # Initial neighbors
  seen_neighbors = set() 
  seen_neighbors.add(node) # Sees itself as seen
  neighbors = set(graph.neighbors(node))
  seen_neighbors.update(neighbors) # Add initial neighbors to seen
  for _ in range(n - 1):
    # New neighbors are the neighbors of the current neighbors
    new_neighbors = set()
    for neighbor in neighbors:
      # Only consider neighbors that haven't been seen yet
      for nn in graph.neighbors(neighbor):
        if nn not in seen_neighbors:
          new_neighbors.add(nn)
    
    # If no unseen neighbors were found, stop
    if not new_neighbors:
      print("No more new neighbors found.")
      break

    # Add to seen neighbors and move frontier forward
    seen_neighbors.update(new_neighbors)
    neighbors = new_neighbors
    
  return neighbors