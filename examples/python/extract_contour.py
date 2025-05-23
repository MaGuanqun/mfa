#### import the simple module from the paraview
from paraview.simple import *
import sys
import time
import numpy as np

file_name = sys.argv[1]
value = float(sys.argv[2])

obj_file_name = sys.argv[3]
point_file_name = sys.argv[4]


print(f"File Name: {file_name}")

start_time_1 = time.time()

# create a new 'Legacy VTK Reader'
s3dmfa_3dvtk = LegacyVTKReader(registrationName='s3d.mfa_3d.vtk', FileNames=[file_name])

tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=s3dmfa_3dvtk)
tetrahedralize1.UpdatePipeline()

start_time_2 = time.time()

contour1 = Contour(registrationName='Contour1', Input=tetrahedralize1)
contour1.ContourBy = ['POINTS', 'var0']
contour1.Isosurfaces = [value]
contour1.PointMergeMethod = 'Uniform Binning'

contour1.UpdatePipeline()

SaveData(obj_file_name, proxy=contour1, PointDataArrays=['var0'])

print(obj_file_name)

end_time = time.time()
execution_time = end_time - start_time_2
total_time = end_time - start_time_1
print(f"Execution Time: {execution_time} seconds")
print(f"Total Time: {total_time} seconds")

# Get the points from the contour
points = contour1.GetClientSideObject().GetOutput().GetPoints().GetData()

# Convert the points to a numpy array
points_array = np.array(points, dtype=np.float64)



# Save the points to a binary file
points_array.tofile(point_file_name)
