from netCDF4 import Dataset
import numpy as np
import gzip
from scipy.ndimage import generic_filter

import meshio


# file_path = '../../build/examples/boussinesq/boussinesq.nc'
# nc_data = Dataset(file_path, mode='r')


# # print(nc_data)

# u0 = nc_data.variables['u'][1999, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1999, :, :].filled(-2.0)

# u1 = nc_data.variables['u'][2000, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v1 = nc_data.variables['v'][2000, :, :].filled(-2.0)



# max_value = np.max(u1)
# min_value = np.min(u1)

# print(max_value)
# print(min_value)


# velocity_t1500 = np.sqrt(u0**2 + v0**2)
# velocity_t1501 = np.sqrt(u1**2 + v1**2)


# # print(velocity_t1500.shape)

# velocity_t1500.astype('float32').tofile('velocity_t2000.bin')
# velocity_t1501.astype('float32').tofile('velocity_t2001.bin')

# nc_data.close()





# # Replace 'your_file.nc' with the path to your NetCDF file
# file_path = '../../build/examples/vortexstreet/cylinder2d.nc'
# nc_data = Dataset(file_path, mode='r')

# u0 = nc_data.variables['u'][1499, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1499, :, :].filled(-2.0)

# u1 = nc_data.variables['u'][1500, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v1 = nc_data.variables['v'][1500, :, :].filled(-2.0)



# max_value = np.max(u1)
# min_value = np.min(u1)

# print(max_value)
# print(min_value)


# velocity_t1500 = np.sqrt(u0**2 + v0**2)
# velocity_t1501 = np.sqrt(u1**2 + v1**2)


# print(velocity_t1500.shape)

# velocity_t1500.astype('float32').tofile('velocity_t1500.bin')
# velocity_t1501.astype('float32').tofile('velocity_t1501.bin')

# nc_data.close()




# amira_mesh = nib.load('../../build/examples/fluid_simu_ML/3000.am')
# data = amira_mesh.get_fdata()
# print(data.shape)
# print(data)




XDIM = 500
YDIM = 500
ZDIM = 100
TDIM = 1  # Single time step for each file

# Load the binary data
file_path = '../../build/examples/hurricane_isabel/Pf30.bin.gz'



with gzip.open(file_path, 'rb') as f:
    # Step 2: Read the data as Big Endian float32
    data = np.frombuffer(f.read(), dtype='>f4')

expected_size = XDIM * YDIM * ZDIM * TDIM
if data.size != expected_size:
    raise ValueError(f"Data size {data.size} does not match expected size {expected_size}.")




data = data.reshape((TDIM,ZDIM,YDIM, XDIM))




slice_z50 = data[0, 50, :, :].copy()  

missing_value = 1e35
slice_z50[slice_z50 == missing_value] = np.nan

def nanmean_filter(values):
    valid_values = values[~np.isnan(values)]
    if valid_values.size > 0:
        return np.mean(valid_values)
    else:
        return np.nan
    
while np.isnan(slice_z50).any():
    slice_z50 = generic_filter(slice_z50, nanmean_filter, size=3, mode='constant', cval=np.nan)

slice_z50 = np.nan_to_num(slice_z50, nan=missing_value)

max_value = np.max(slice_z50)
min_value = np.min(slice_z50)

print(max_value)
print(min_value)


# print(np.ma.isMaskedArray(slice_z50))

slice_z50.astype('float32').tofile('../../build/examples/hurricane_isabel/Pf30_50.bin')


print("Shape of slice_z50:", slice_z50.shape)
max_value = np.max(slice_z50)
min_value = np.min(slice_z50)

print(max_value)
print(min_value)

has_nan = np.isnan(slice_z50).any()

# Print the result
if has_nan:
    print("There are NaN values in slice_z50.")
else:
    print("There are no NaN values in slice_z50.")
