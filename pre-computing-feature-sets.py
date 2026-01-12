from synkit.IO import rsmi_to_graph
from _blake2 import blake2b
from networkx.algorithms import all_pairs_shortest_path
from drf_implementation import vertex_drf, edge_drf, shortest_path_drf
import pandas as pd

def calculate_symmetric_difference_off_dict(dict1, dict2):
  # Calculate for each key in both dicts the symmetric difference of their value sets
  symm_diff = dict()
  all_keys = set(dict1.keys()).union(set(dict2.keys()))
  for key in all_keys:
    set1 = set(dict1.get(key, []))
    set2 = set(dict2.get(key, []))
    symm_diff[key] = set1.symmetric_difference(set2)
  return symm_diff

# Load schneider50k-clean.tsv to dataframe
df = pd.read_csv("schneider50k_clean.tsv",sep="\t")

# For each reaction in the dataframe, compute the phi_vertex_dict_graph, phi_edge_graph, and phi_shortest_path_graph
# Store in a new dataframe consisting of 
### educt_phi_vertex_dict, product_phi_vertex_dict, symmetric_difference_vertex_dict,
### educt_phi_edge, product_phi_edge, symmetric_difference_edge,
### educt_phi_shortest_path, product_phi_shortest_path, symmetric_difference_shortest_path
df_features = pd.DataFrame(columns=[
    "educt_phi_vertex_dict", "product_phi_vertex_dict", "symmetric_difference_vertex_dict",
    "educt_phi_edge", "product_phi_edge", "symmetric_difference_edge",
    "educt_phi_shortest_path", "product_phi_shortest_path", "symmetric_difference_shortest_path"
])

for index, row in df.iterrows():
  rsmi = row["clean_rxn"]

  educt_phi_vertex_dict, product_phi_vertex_dict, symm_diff_vertices_test = vertex_drf(rsmi)

  educt_phi_edge, product_phi_edge, symm_diff_edges_test = edge_drf(rsmi)

  educt_phi_shortest_path, product_phi_shortest_path, symm_diff_shortest_paths_test = shortest_path_drf(rsmi)

  # Print short summary
  print(f"Processed reaction {index}/{len(df)}")
  df_features = pd.concat([df_features, pd.DataFrame([{
      "educt_phi_vertex_dict": educt_phi_vertex_dict,
      "product_phi_vertex_dict": product_phi_vertex_dict,
      "symmetric_difference_vertex_dict": symm_diff_vertices_test,
      "educt_phi_edge": educt_phi_edge,
      "product_phi_edge": product_phi_edge,
      "symmetric_difference_edge": symm_diff_edges_test,
      "educt_phi_shortest_path": educt_phi_shortest_path,
      "product_phi_shortest_path": product_phi_shortest_path,
      "symmetric_difference_shortest_path": symm_diff_shortest_paths_test
  }])], ignore_index=True)

  # Test: Break after 1 reaction
  if index >= 0:
    # Write to excel
    #df_features.to_excel("schneider50k_features.xlsx", index=False)
    break;