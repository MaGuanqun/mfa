import pyvista as pv
import numpy as np
import argparse
#use "ttk" as environment

def add_z(file_path, file_name):
    # mesh_path = '../../build/examples/cesm/up_50.ply'     # Update this path
    mesh = pv.read(file_path + file_name)
    points = mesh.points
    
    max_z=np.max(points[:,2])+5
    
    min_x = np.min(points[:,0])
    max_x = np.max(points[:,0])
    min_y = np.min(points[:,1])
    max_y = np.max(points[:,1])
    
    vertices = np.array([[min_x, min_y, max_z],[min_x, max_y, max_z],[max_x, max_y, max_z],[max_x, min_y, max_z]])
    faces = np.hstack([[4, 0, 1, 2, 3]]) 
    square_mesh = pv.PolyData(vertices, faces)
    
    # mesh.points = points    
    
    square_mesh.save(file_path +'add_z_'+file_name)



parser = argparse.ArgumentParser(description='add z value.')
parser.add_argument('-n', '--name', type=str, default='cesm', help='object name')
parser.add_argument('-p', '--path', type=str, default='ces,', help='folder path to all files')

args = parser.parse_args()

obj_name=args.name
file_path = args.path

k=[10,11]
for i in k:
    add_z(file_path, obj_name+'_'+str(i)+'.ply')
    


# # Criteria for filtering: vertices within a certain range
# # Example: x in [xmin, xmax], y in [ymin, ymax], z in [zmin, zmax]
# xmin, xmax = 1229.3, 1229.5
# ymin, ymax = 932, 932.2

# condition = (mesh.points[:, 0] >= xmin) & (mesh.points[:, 0] <= xmax) & (mesh.points[:, 1] >= ymin) & (mesh.points[:, 1] <= ymax) 
# # Find vertices within the specified range
# submesh = mesh.extract_points(condition)
# # submesh = submesh.extract_points(submesh.points[:, 0] <= xmax)
# # submesh = submesh.extract_points(submesh.points[:, 1] >= ymin)
# # submesh = submesh.extract_points(submesh.points[:, 1] <= ymax)


# # Save or visualize the extracted submesh
# submesh.save('../../build/examples/cesm/extracted_submesh.vtk')  # Save the submesh to a file
# submesh.plot(show_edges=True)  # Visualize the submesh