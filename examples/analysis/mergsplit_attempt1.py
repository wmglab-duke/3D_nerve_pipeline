import matplotlib.pyplot as plt
import networkx as nx

# Expanded synthetic data with more fascicles and complex merges and splits
edges = [
    (0, 1, 1),
    (0, 2, 1),  # Split from 0 to 1 and 2
    (1, 3, 2),
    (2, 3, 2),  # Merge from 1 and 2 to 3
    (3, 4, 3),
    (3, 5, 3),  # Split from 3 to 4 and 5
    (4, 6, 4),
    (5, 6, 4),  # Merge from 4 and 5 to 6
    (1, 7, 2),
    (7, 8, 3),
    (2, 8, 3),  # More complex branching and merging
    (6, 9, 5),
    (8, 9, 5),  # Merge into 9
    (9, 10, 6),
    (9, 11, 6),  # Split from 9 to 10 and 11
    (10, 12, 7),
    (11, 13, 7),  # Split further
    (4, 14, 4),
    (14, 15, 5),  # Independent branch out from 4
    (15, 12, 6),  # Joining into 12
    (12, 16, 8),
    (13, 16, 8),  # Final merge
    (5, 17, 4),
    (17, 18, 5),
    (11, 18, 6),  # Additional branches and merge
]

# Create a directed graph
G = nx.DiGraph()
positions = {}
levels = {}

for src, tgt, level in edges:
    G.add_edge(src, tgt)
    levels[src] = level - 1
    levels[tgt] = level

# Assigning positions for clear visualization of merges and splits
max_width = max(levels.values())
positions = {node: (levels[node], 0) for node in G.nodes()}
for level in range(max_width + 1):
    nodes_at_level = [node for node in positions if positions[node][0] == level]
    for i, node in enumerate(sorted(nodes_at_level)):
        positions[node] = (level, i - len(nodes_at_level) / 2)

# Drawing the graph
plt.figure(figsize=(14, 10))
nx.draw(
    G,
    pos=positions,
    node_size=500,
    node_color='skyblue',
    with_labels=True,
    font_weight='bold',
    arrowstyle='-|>',
    arrowsize=20,
)
plt.title('Complex Nerve Fascicle Topology Showing Merges and Splits')
plt.gca().invert_yaxis()  # Invert y-axis to show the root at the top
plt.show()
