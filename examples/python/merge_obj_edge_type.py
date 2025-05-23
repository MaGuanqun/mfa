from paraview.simple import *
import pandas as pd
import os
import csv
import argparse
import time

def read_obj(file_path):
    vertices = []
    edges = []
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split()
            if parts[0] == "v":  # Vertex
                vertices.append(tuple(map(float, parts[1:])))
            elif parts[0] == "l":  # Edge
                edges.append(tuple(map(lambda x: int(x) - 1, parts[1:])))  # Convert to 0-based indexing
    return vertices, edges

def read_edge_values(file_path):
    with open(file_path, "r") as file:
        return [int(line.strip()) for line in file]
    
def create_programmable_source(vertices, edges, edge_values):
    """Create a Programmable Source in ParaView."""
    programmableSource = ProgrammableSource()

    # Set the dataset type
    programmableSource.OutputDataSetType = "vtkPolyData"

    # Generate the script for the Programmable Source
    vertices_str = str(vertices)
    edges_str = str(edges)
    edge_values_str = str(edge_values)

    programmableSource.Script = f"""
# Define vertices, edges, and edge values
vertices = {vertices_str}
edges = {edges_str}
edge_values = {edge_values_str}

# Import necessary modules from ParaView
from paraview.vtk import vtkPoints, vtkCellArray, vtkFloatArray, vtkLine

# Create points
points = vtkPoints()
for vertex in vertices:
    points.InsertNextPoint(vertex)

# Create lines
lines = vtkCellArray()
scalars = vtkFloatArray()
scalars.SetName("EdgeValues")
for i, edge in enumerate(edges):
    line = vtkLine()
    line.GetPointIds().SetId(0, edge[0])
    line.GetPointIds().SetId(1, edge[1])
    lines.InsertNextCell(line)
    scalars.InsertNextValue(edge_values[i])

# Create the output polydata
output.SetPoints(points)
output.SetLines(lines)
output.GetCellData().SetScalars(scalars)
"""
    programmableSource.UpdatePipeline()
    return programmableSource


def main(obj_file, edge_values_file, output_file):
    # Read the input files
    vertices, edges = read_obj(obj_file)
    edge_values = read_edge_values(edge_values_file)

    assert len(edges) == len(edge_values), "Mismatch between edge count and edge values!"

    # Create the dataset in ParaView
    programmableSource = create_programmable_source(vertices, edges, edge_values)

    # Save the output as a .vtp file
    SaveData(output_file, proxy=programmableSource)
    print(f"File saved: {output_file}")


parser = argparse.ArgumentParser(description='merge obj edge types.')


parser.add_argument('-i', '--input_obj_name', type=str, default='file_name.obj', help='obj file')
parser.add_argument('-j', '--input_edge_type', type=str, default='file_name.txt', help='edge value file')
parser.add_argument('-o','--output_name', type=str, default='output.vtp', help='output csv name')


args = parser.parse_args()

obj_file_name=args.input_obj_name
type_file_name=args.input_edge_type
output_file=args.output_name

main(obj_file_name, type_file_name, output_file)


