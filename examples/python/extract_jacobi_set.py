#### import the simple module from the paraview
from paraview.simple import *
import sys
import time
import numpy as np
import os

def load_obj(filename):
    vertices = []
    edge_lines = []
    header_lines = []
    with open(filename, 'r') as f:
        for line in f:
            # Save header or other lines if needed
            if line.startswith("v "):
                parts = line.strip().split()
                # Parse x,y,z coordinates
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith("l ") or (line.startswith("f ") and len(line.strip().split())==3):
                # Assume "l" or "f" lines with exactly 2 vertex indices define an edge.
                edge_lines.append(line)
            else:
                header_lines.append(line)
    vertices = np.array(vertices)
    return vertices, edge_lines, header_lines

def merge_duplicate_vertices(vertices, tol=1e-4):
    """
    Merges duplicate vertices using a tolerance.
    Returns an array of unique vertices and a mapping from original index to new index.
    """
    # Multiply coordinates by factor = 1/tol and round to integers
    factor = 1.0 / tol
    rounded = np.round(vertices * factor).astype(np.int64)
    # Use np.unique along axis 0.
    unique_rounded, unique_indices, inverse_indices = np.unique(rounded, axis=0, return_index=True, return_inverse=True)
    unique_vertices = vertices[unique_indices]
    return unique_vertices, inverse_indices


def update_edge_indices(edge_lines, mapping):
    """
    Updates edge definitions to use new vertex indices.
    The 'mapping' array maps old indices to new indices.
    OBJ vertex indices are 1-indexed.
    """
    updated_edges = []
    for line in edge_lines:
        parts = line.split()
        try:
            # Convert OBJ indices (1-indexed) to 0-indexed
            idx0 = int(parts[1]) - 1
            idx1 = int(parts[2]) - 1
        except (IndexError, ValueError):
            continue
        # Map to the new indices and convert back to 1-indexed.
        new_idx0 = mapping[idx0] + 1
        new_idx1 = mapping[idx1] + 1
        updated_edges.append(f"{parts[0]} {new_idx0} {new_idx1}")
    return updated_edges

def compute_bounds(vertices):
    xmin = vertices[:,0].min()
    xmax = vertices[:,0].max()
    ymin = vertices[:,1].min()
    ymax = vertices[:,1].max()
    zmin = vertices[:,2].min()
    zmax = vertices[:,2].max()
    return xmin, xmax, ymin, ymax, zmin, zmax

def is_on_boundary(vertex, bounds, tol=1e-4):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    x, y, z = vertex
    return (abs(x - xmin) < tol or abs(x - xmax) < tol or
            abs(y - ymin) < tol or abs(y - ymax) < tol)

def filter_edges(vertices, edge_lines, tol=1e-6):
    """
    Returns the kept edge lines and a set of vertex indices used by kept edges.
    """
    bounds = compute_bounds(vertices)
    kept_edges = []
    used_indices = set()
    for line in edge_lines:
        parts = line.split()
        try:
            # OBJ vertex indices are 1-indexed; convert to 0-index.
            idx0 = int(parts[1]) - 1
            idx1 = int(parts[2]) - 1
        except (IndexError, ValueError):
            # Skip lines that don't conform.
            continue
        v0 = vertices[idx0]
        v1 = vertices[idx1]
        # Check if both endpoints lie on the boundary.
        if is_on_boundary(v0, bounds, tol) and is_on_boundary(v1, bounds, tol):
            # Skip this edge.
            continue
        else:
            kept_edges.append((idx0, idx1, parts[0]))  # keep the type ("l" or "f")
            used_indices.add(idx0)
            used_indices.add(idx1)
    return kept_edges, used_indices

def save_obj2(filename, vertices, edge_lines, header_lines):
    with open(filename, 'w') as f:
        # Write any header or comment lines.
        for line in header_lines:
            f.write(line + "\n")
        # Write vertex definitions.
        for v in vertices:
            f.write("v {:.6f} {:.6f} {:.6f}\n".format(v[0], v[1], v[2]))
        # Write edge definitions.
        for line in edge_lines:
            f.write(line + "\n")

def save_obj(filename, vertices, kept_edges, header_lines,used_indices):
    """
    Saves an OBJ file with only the vertices that are used by kept_edges.
    The vertices are re-indexed, and the edge definitions are updated.
    """
    
    # Build a sorted list of used vertex indices.
    used_list = sorted(list(used_indices))
    
    # print(used_list)
    
    # Create a mapping from old index to new index.
    mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(used_list)}
    # Create a new list of vertices (only the used ones).
    new_vertices = vertices[used_list]

    with open(filename, 'w') as f:
        # Write header lines first if any.
        for line in header_lines:
            f.write(line)
        # Write new vertex definitions.
        for v in new_vertices:
            f.write("v {:.6f} {:.6f} {:.6f}\n".format(v[0], v[1], v[2]))
        # Write the kept edges with re-indexed vertex indices.
        for edge in kept_edges:
            idx0, idx1, edge_type = edge
            new_idx0 = mapping[idx0] + 1  # convert back to 1-indexed
            new_idx1 = mapping[idx1] + 1
            # Write using the original edge keyword (e.g., "l" or "f")
            f.write(f"{edge_type} {new_idx0} {new_idx1}\n")
    return new_vertices




file_name = sys.argv[1]

obj_file_name = sys.argv[2]
point_file_name = sys.argv[3]


print(f"File Name: {file_name}")

start_time_1 = time.time()

# create a new 'Legacy VTK Reader'
s3dmfa_3dvtk = LegacyVTKReader(registrationName='s3d.mfa_3d.vtk', FileNames=[file_name])

tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=s3dmfa_3dvtk)
tetrahedralize1.UpdatePipeline()


# Free the memory of s3dmfa_3dvtk
del s3dmfa_3dvtk

# SaveData(file_name+".vtk", proxy=tetrahedralize1)

start_time_2 = time.time()

# contour1 = Contour(registrationName='Contour1', Input=tetrahedralize1)
# contour1.ContourBy = ['POINTS', 'var0']
# contour1.Isosurfaces = [value]
# contour1.PointMergeMethod = 'Uniform Binning'


# create a new 'TTK JacobiSet'
jacobi_set = TTKJacobiSet(registrationName='TTKJacobiSet1', Input=tetrahedralize1)
jacobi_set.UComponent = ['POINTS', 'var0']
jacobi_set.VComponent = ['POINTS', 'var1']
jacobi_set.UOffsetField = ['POINTS', 'var0']
jacobi_set.VOffsetField = ['POINTS', 'var1']


jacobi_set.UpdatePipeline()



# Save the resulting data (now without edges having both endpoints on the boundary)
SaveData(obj_file_name, proxy=jacobi_set)

end_time = time.time()

obj_data = open(obj_file_name, 'r').read()
if os.path.exists(obj_file_name):
    os.remove(obj_file_name)

# Substitute 'f' with 'e'
obj_data = obj_data.replace('f', 'l')

# Save the modified OBJ data to a file
with open(obj_file_name, 'w') as obj_file:
    obj_file.write(obj_data)



vertices, edge_lines, header_lines = load_obj(obj_file_name)
# Delete the obj_file_name file
if os.path.exists(obj_file_name):
    os.remove(obj_file_name)
kept_edges, used_indices  = filter_edges(vertices, edge_lines, tol=1e-3)
new_vertices= save_obj(obj_file_name, vertices, kept_edges, header_lines,used_indices)


vertices, edge_lines, header_lines = load_obj(obj_file_name)
# Delete the obj_file_name file
if os.path.exists(obj_file_name):
    os.remove(obj_file_name)
unique_vertices, mapping = merge_duplicate_vertices(vertices, tol=1e-4)
updated_edges = update_edge_indices(edge_lines, mapping)

save_obj2(obj_file_name, unique_vertices, updated_edges, header_lines)
# Read the OBJ file

# Print the modified OBJ data

print(obj_file_name)

execution_time = end_time - start_time_2
total_time = end_time - start_time_1
print(f"Execution Time: {execution_time} seconds")
print(f"Total Time: {total_time} seconds")


# vtk_output = servermanager.Fetch(jacobi_set)

# # Extract the points data
# points = vtk_output.GetPoints().GetData()



# Convert the points to a numpy array
points_array = np.array(new_vertices, dtype=np.float64)

print(points_array)


# Save the points to a binary file
points_array.tofile(point_file_name)
