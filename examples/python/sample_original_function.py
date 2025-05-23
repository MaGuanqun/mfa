import numpy as np

# Define the 1D function with special handling for t == 0
def f(t):
    tol = 1e-8
    return np.where(np.abs(t) > tol, np.sin(5*t)/t, 5.0)


def F(x, y, s=1.0):
    a = 418.9829
    d = 2  # dimension of the domain (x and y)
    # Compute the term: x*sin(sqrt(|x|)) + y*sin(sqrt(|y|))
    term = x * np.sin(np.sqrt(np.abs(x))) + y * np.sin(np.sqrt(np.abs(y)))
    return s * 0.5 * (a * d - term)


# Set up the grid parameters
nx, ny = 301, 301  # grid resolution in x and y directions
# x = np.linspace(-10.5* 10.5 * np.pi* np.pi, 10.5* 10.5 * np.pi* np.pi, nx)
# y = np.linspace(-10.5* 10.5 * np.pi* np.pi, 10.5* 10.5 * np.pi* np.pi, ny)

x = np.linspace(-2.0 * np.pi, 2.0 * np.pi, nx)
y = np.linspace(-2.0 * np.pi, 2.0 * np.pi, ny)

# Create a 2D grid. Using 'xy' indexing so that X.shape and Y.shape are (ny, nx)
X, Y = np.meshgrid(x, y, indexing='xy')
# Compute the function value at each grid point: F(x,y) = f(x) + f(y)
Z = f(X) + f(Y)

# Z= f(X, Y, s=1.0)


# Build the list of vertices. Each vertex will be (x, y, z)
# Note: the grid is traversed row-by-row (y index first)
vertices = []
for j in range(ny):
    for i in range(nx):
        vertices.append([x[i], y[j], Z[j, i]])
vertices = np.array(vertices)

# Build triangle connectivity by splitting each grid cell into 2 triangles.
# For a cell defined by indices:
#   v0 = (i,   j)
#   v1 = (i+1, j)
#   v2 = (i,   j+1)
#   v3 = (i+1, j+1)
# We create two triangles: (v0, v1, v3) and (v0, v3, v2).
triangles = []
for j in range(ny - 1):
    for i in range(nx - 1):
        # Compute the vertex indices in the flattened vertices array:
        # The index ordering is: index = i + j * nx.
        v0 = i + j * nx
        v1 = (i + 1) + j * nx
        v2 = i + (j + 1) * nx
        v3 = (i + 1) + (j + 1) * nx
        # First triangle: vertices v0, v1, v3
        triangles.append([v0, v1, v3])
        # Second triangle: vertices v0, v3, v2
        triangles.append([v0, v3, v2])

# Write the data to a legacy VTK PolyData file with triangles
vtk_filename = "sinc_2d_triangles.vtk"
with open(vtk_filename, "w") as f_out:
    f_out.write("# vtk DataFile Version 3.0\n")
    f_out.write("Triangulated sinc_sum 2D data\n")
    f_out.write("ASCII\n")
    f_out.write("DATASET POLYDATA\n")
    # Write the points
    f_out.write("POINTS {} float\n".format(len(vertices)))
    for pt in vertices:
        f_out.write("{} {} {}\n".format(pt[0], pt[1], pt[2]))
    # Write the triangle connectivity.
    # For POLYGONS, each triangle is stored with 4 numbers: the count (3) and then 3 point indices.
    num_triangles = len(triangles)
    total_size = num_triangles * 4
    f_out.write("POLYGONS {} {}\n".format(num_triangles, total_size))
    for tri in triangles:
        f_out.write("3 {} {} {}\n".format(tri[0], tri[1], tri[2]))
    # Add the scalar function values as point data.
    f_out.write("POINT_DATA {}\n".format(len(vertices)))
    f_out.write("SCALARS function float 1\n")
    f_out.write("LOOKUP_TABLE default\n")
    # The function values are written in the same order as the vertices.
    for value in Z.flatten():
        f_out.write("{}\n".format(value))

print("VTK file '{}' has been written.".format(vtk_filename))
