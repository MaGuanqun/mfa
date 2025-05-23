import pyvista as pv
import numpy as np

# min=[48,48,46]
# max=[68,68,80]
# min=[24,24,68]
# max=[44,44,103]
# min=[41,7,34]
# max=[61,27,68]




def save3D(min,max,name,index):
    points = np.array([
        [min[0], min[1], min[2]], # Point 0
        [max[0], min[1], min[2]], # Point 1
        [max[0], max[1], min[2]], # Point 2
        [min[0], max[1], min[2]], # Point 3
        [min[0], min[1], max[2]], # Point 4
        [max[0], min[1], max[2]], # Point 5
        [max[0], max[1], max[2]], # Point 6
        [min[0], max[1], max[2]]  # Point 7
    ], dtype=float)

    faces = np.array([
        np.hstack([[4, 0, 3, 7, 4]]),  # Back face
        np.hstack([[4, 0, 4, 5, 1]]),  # Bottom face
        np.hstack([[4, 3, 2, 6, 7]]),  # Left face
        np.hstack([[4, 1, 5, 6, 2]]),  # Front face
        np.hstack([[4, 4, 7, 6, 5]]),  # Top face
        np.hstack([[4, 0, 1, 2, 3]])   # Right face
    ]).flatten()
    cube = pv.PolyData(points, faces)
    # Save to a PLY file
    filename = "../build/examples/"+name+"/"+name+"_"+index+".ply"
    cube.save(filename)

def save2D(min,max,name,index):
    points = np.array([[min[0], min[1], 2], # Point 0
        [max[0], min[1], 2], # Point 1
        [max[0], max[1], 2], # Point 2
        [min[0], max[1], 2] # Point 3
        ], dtype=float)
    faces = np.array([
        np.hstack([[4, 0, 1, 2, 3]])   # Right face
    ]).flatten()
    square = pv.PolyData(points, faces)
    
    # Save to a PLY file
    filename = "../build/examples/"+name+"/"+name+"_"+index+"_block.ply"
    
    print(filename)
    
    square.save(filename)
    
min=[]
max=[]


# min=[48,48,46]
# max=[68,68,80]
# min=[24,24,68]
# max=[44,44,103]
# min.append([41,7,34])
# max.append([61,27,68])


# min.append([51.9, 48, 51])
# max.append([54.9, 51, 52.6])

# min.append([29.8, 29.7, 88.7])
# max.append([32.9, 32.7, 90.3])

# min.append([58, 11.5, 58.1])
# max.append([61, 14.5, 59.7])

# min.append([74.65, 75.74,  94.1])
# max.append([76.08, 78.29, 96.65])

# min.append([196.84, 196.74 ])
# max.append([217.93, 212.91])

# min.append([158.18,  99.72])
# max.append([179.27, 115.89])

# min.append([45.62,  60.69, 206.55])
# max.append([ 47.05, 63.24, 209.1])


min.append([-0.207732, 0.350515 ])
max.append([-0.175773, 0.382474])

name="expotential"

index=[0]

for i in index:
    save2D(min[i],max[i],name,str(i+5))