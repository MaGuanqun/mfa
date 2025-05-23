from paraview.simple import *

import sys

file_name = sys.argv[1]

print(f"File Name: {file_name}")

# create a new 'Wavefront OBJ Reader'
obj = WavefrontOBJReader(registrationName='isocontour.obj', FileName=file_name)


connectivity = Connectivity(registrationName='Connectivity1', Input=obj)

connectivity.UpdatePipeline()


# Extract RegionId array from the output of the Connectivity filter
connectivity_data = servermanager.Fetch(connectivity)  # Fetch data to Python
region_array = connectivity_data.GetPointData().GetArray("RegionId")

# Get unique region IDs (number of connected components)
num_components = region_array.GetRange()[1] + 1  # Max RegionId + 1 gives the count


# Get vertex and edge counts
info = obj.GetDataInformation()
num_vertices = info.GetNumberOfPoints()
num_edges = info.GetNumberOfCells()  # Since cells are edges


data = servermanager.Fetch(obj)
# Create a dictionary to store node degrees
node_degrees = {i: 0 for i in range(num_vertices)}


# Iterate over edges (cells) and count connections
for i in range(num_edges):
    edge = data.GetCell(i)  # Get the edge
    point_ids = [edge.GetPointId(j) for j in range(edge.GetNumberOfPoints())]
    
    # Increment the degree count for each connected vertex
    for pid in point_ids:
        node_degrees[pid] += 1


num_single_points = sum(1 for degree in node_degrees.values() if degree == 0)
print(f"Single Points: {num_single_points}")




# Compute the number of loops (cycles)
num_loops = num_edges - num_vertices + num_components + num_single_points

# Print results
print(f"Connected Components: {num_components}")
print(f"Loops (Cycles): {num_loops}")





# Find nodes with 3 or more edges
high_degree_nodes = [pid for pid, degree in node_degrees.items() if degree >= 3]
print("nodes with 3+ edges:", len(high_degree_nodes))

# points = data.GetPoints()
# for pid in high_degree_nodes:
#     print(points.GetPoint(pid))

# print(high_degree_nodes)