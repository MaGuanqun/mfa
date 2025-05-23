import numpy as np
from scipy.ndimage import zoom

x,y,z=288,512,512
num_elements=x*y*z

dtype=np.float32

file = "../../build/examples/ori_data/dd07g_xxsmall_le.xyz"


output_file = "../../build/examples/ori_data/rti_144_256_256.xyz"

with open(file, 'rb') as file:
    data = np.fromfile(file, dtype=dtype, count=num_elements*3)

data_shaped = data.reshape(-1,3)
norms=np.linalg.norm(data_shaped, axis=1)

print(norms)

data = norms.reshape((x, y, z), order='F') #F  For column-major order
                #('C')  # For row-major order (default)

# print(data)         
                
# data_block = data[0:72, 0:80, 0:80]


fx, fy, fz = 0.5, 0.5, 0.5
downsampled_data = zoom(data, (fx, fy, fz))
print(downsampled_data.shape)


downsampled_data = downsampled_data.flatten('F').astype(np.float32)
# downsampled_data_fortran = np.asfortranarray(data_block)
print(downsampled_data)
print(downsampled_data.shape)

with open(output_file, 'wb') as file:
    downsampled_data.tofile(file)
    