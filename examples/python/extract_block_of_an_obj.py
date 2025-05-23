import sys


input_file = sys.argv[1]
output_file = sys.argv[2]
block_str = sys.argv[3]

def is_in_block(x, y, block):
    """
    Check if the point (x,y,z) is within the block.
    Block is defined as a tuple: (xmin, xmax, ymin, ymax, zmin, zmax)
    """
    xmin, xmax, ymin, ymax = block
    return (xmin <= x <= xmax and ymin <= y <= ymax)



block = block_str.split("-")
block = list(map(float, block))


# Read vertices and edges from the input file.
vertices = []  # list of tuples (x, y, z)
edges = []     # list of lists of vertex indices (1-indexed as in OBJ files)
with open(input_file, 'r') as f:
    for line in f:
        # Process vertex lines
        if line.startswith("v "):
            parts = line.strip().split()
            # Convert coordinate strings to floats
            x, y, z = map(float, parts[1:4])
            vertices.append((x, y, z))
        # Process edge lines (starting with 'l')
        elif line.startswith("l "):
            parts = line.strip().split()
            # Convert indices to integers (OBJ indices are 1-based)
            indices = list(map(int, parts[1:]))
            edges.append(indices)

# Filter out vertices that are in the specified block.
# Build a mapping from the original vertex index to a new index.
new_vertices = []
old_to_new = {}  # key: old index, value: new index
for i, (x, y, z) in enumerate(vertices, start=1):
    if is_in_block(x, y, block):
        # Only keep vertices not in the block.
        old_to_new[i] = len(new_vertices) + 1  # new index (1-based)
        new_vertices.append((x, y, z))

# Filter edges:
# Keep an edge only if all its vertex indices are present in the new vertex mapping.
new_edges = []
for edge in edges:
    new_edge = []
    valid_edge = True
    for idx in edge:
        if idx in old_to_new:
            new_edge.append(old_to_new[idx])
        else:
            valid_edge = False
            break  # one vertex is removed; skip the whole edge
    if valid_edge:
        new_edges.append(new_edge)

# Write the filtered vertices and edges to a new OBJ file.
with open(output_file, 'w') as f:
    # Write vertices
    for (x, y, z) in new_vertices:
        f.write(f"v {x} {y} {z}\n")
    # Write edges
    for edge in new_edges:
        edge_str = " ".join(map(str, edge))
        f.write(f"l {edge_str}\n")

print(f"Filtered OBJ file saved as {output_file}")
