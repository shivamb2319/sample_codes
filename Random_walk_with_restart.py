#importing the necessary libraries

import numpy as np
import pandas as pd
import networkx as nx

# reading the interactome file

df= pd.read_excel('your_file.xlsx') 

# removing the self loops 

df = df[df['Gene_A'] != df['Gene_B']]

# creating a graph by dropping a edge from node in Gene_A column to node in Gene_B column

G = nx.from_pandas_edgelist(df, 'Gene_A', 'Gene_B')

# extracting the largest connected component from the formed graph

largest_component = max(nx.connected_components(G), key=len)

# updating the dataframe (df) to contain only the nodes from the largest connected component

df = df[df.apply(lambda row: row['Gene_A'] in largest_component and row['Gene_B'] in largest_component, axis=1)]

# creating an adjacency matrix from the filtered dataframe, 1 for every connection between 2 nodes and 0 for no connection

num_nodes1=list(df['Gene_A'])
num_nodes2=list(df['Gene_B'])
sum=num_nodes1+num_nodes2

unique_list = []
for item in sum:
    if item not in unique_list:
        unique_list.append(item)
A=unique_list

num_nodes = len(A)
adj_matrix = pd.DataFrame(0, index=A, columns=A)

for row in df.itertuples(index=False):
    adj_matrix.at[row.Gene_A, row.Gene_B] = 1
    adj_matrix.at[row.Gene_B, row.Gene_A] = 1


# creating a diagonal matrix from the adjacency matrix and then normalized adjacency matrix

col_sums = np.sum(adj_matrix, axis=0)

reciprocal_sqrt_col_sums = np.where(col_sums != 0, 1.0 / np.sqrt(col_sums), 0)

D_matrix = np.diag(reciprocal_sqrt_col_sums)

n_adj_matrix= (D_matrix)@(adj_matrix)@(D_matrix)

# sum2one fnc for normalizing the matrix ( sum of the column elements should be 1) 

def sum2one(x):
     global p0_matrix2
     column_sums = np.sum(x, axis=0)
     p0_matrix2 = (x)*(1/ column_sums)
     return p0_matrix2

def sum2one2(x):
     global normalized_matrix
     column_sums = np.sum(x, axis=0)
     normalized_matrix = x / column_sums
     sum_of_elements = np.sum(normalized_matrix, axis=0)
     adjustment = 1.0 - sum_of_elements
     num_columns = x.shape[1]
     normalized_matrix += adjustment / num_columns
     return normalized_matrix


# RWR code

def random_walk_with_restart(n_adj_matrix, num_steps, restart_prob,epsilon):   
    pt = p0_matrix2.copy()
    global affinity_scores
    affinity_scores = [pt]
    
    for t in range(1, num_steps + 1):
        pt_new = (1 - restart_prob) * np.dot(n_adj_matrix, pt) + restart_prob * p0_matrix2
        np.set_printoptions(precision=8)
        # print(f"The affinity score for step {t} is: {pt_new}")
        affinity_scores.append(pt_new)
        if np.all(np.abs(pt_new - pt) < epsilon):
            break
        pt = pt_new

disease_data = {}
with open("your_file_with_seeded_genes.txt" , "r") as disease_file:
    for line in disease_file:
        line = line.strip()
        if line:
            parts = line.split(":")
            disease_name = parts[0].strip()
            seeded_nodes = [int(node.strip()) for node in parts[1].split(",")]
            disease_data[disease_name] = seeded_nodes

top_nodes_per_disease = {}

for disease_name, seeded_nodes in disease_data.items():
    p0_matrix1 = np.zeros((len(A), 1), dtype=int)
    for i, node in enumerate(A):
        if node in seeded_nodes:
            p0_matrix1[i, 0] = 1

    sum2one(p0_matrix1)
    random_walk_with_restart(n_adj_matrix, 50, 0.75, 1e-6)
    
    pt_matrix = affinity_scores[len(affinity_scores) - 1]
    sum2one2(pt_matrix)

    # Sort the nodes by score in descending order
    sorted_nodes = normalized_matrix[:, 0].argsort()[::-1]
    top_160_nodes = sorted_nodes[:160]

    # Store the top 160 nodes for the current disease
    top_nodes_per_disease[disease_name] = top_160_nodes

with open("output.txt", "w") as result_file:
    for disease_name, top_nodes in top_nodes_per_disease.items():
        result_file.write(f"{disease_name}: {', '.join(str(A[node]) for node in top_nodes)}\n")

print("Top 160 nodes for each disease saved to 'output.txt'")


