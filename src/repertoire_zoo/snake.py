"""
Module Name: Snake

**S**ubject **N**etwork **A**nalysis of **K**ey **E**vents
TODO: Change name to SPIDER to represent the web-like characteristics of a graph network.
Subject Pathway Inference through Dynamic Event Representation?
Subject Pattern Indeicication for Data-drive Event Relations?

Contains functions that create NetworkX graphs.
"""
import airr
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

def create_time_point_graph(
        repertoires,
        subj_IDs,
        pos=None,
        title='Connections Between Time Points'
    ):
    '''
    Parameters
    ----------
    repertoires : dict
        - give this data: airr.read_airr('repertoires.airr.json')
    subj_IDs : list
        - A list of subject IDs to plot between
    pos : dict
        - A dictionary for node positions (Default: None. Results in nx.circular_layout(G))
    title : str
        - Title for the graph (Default: Connections Between Timepoints)
    '''
    df = pd.DataFrame()

    for rep in repertoires['Repertoire']:
        if rep['subject']['subject_id'] in subj_IDs:
            new_df = pd.DataFrame({
                'id' : [rep['subject']['subject_id']], 
                'timepoint' : [rep['sample'][0]['collection_time_point_reference']]
            })
            df = pd.concat([df, new_df], ignore_index=True)

    df = df.sort_values(by=['id', 'timepoint'], ignore_index=True)

    # ID with unique timepoint following current encounter
    connections = {}
    for i in range(len(df)-1):
        if df.loc[i, 'id'] != df.loc[i+1, 'id']:
            pass
        else:
            key = (df.loc[i, 'timepoint'], df.loc[i+1, 'timepoint'])
            if key not in connections:
                connections[key] = 1
            else:
                connections[key] += 1

    # Create graph visualization
    G = nx.Graph()
    times = df['timepoint'].unique()
    G.add_nodes_from(times)

    # add edges
    for (tp1, tp2), w in connections.items():
        G.add_edge(tp1, tp2, weight=w)

    # create layout if not defined
    if pos is None:
        # pos = {tp: (i, 0) for i, tp in enumerate(times)} # linear layout
        # pos = nx.spring_layout(G, seed=42)
        pos = nx.circular_layout(G)

    # node weight
    node_counts = df.groupby('timepoint')['id'].nunique()
    for node in G.nodes():
        if G.has_edge(node, node):
            node_counts[node] += G[node][node]['weight']
    node_sizes = [node_counts.get(tp, 1) * 100 for tp in G.nodes()]
    node_labels = {tp : f"{tp}\n(n={node_counts.get(tp, 1)})" for tp in G.nodes()}

    # edge weight (use widths for normalization for max width)
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    edge_labels = {(u, v): G[u][v]['weight'] for u, v in G.edges()}

    plt.figure(figsize=(10, 6))
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, alpha=0.3)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8)
    nx.draw_networkx_edges(G, pos, width=weights, alpha=0.3)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title(title)
    plt.axis('off')
    plt.show()