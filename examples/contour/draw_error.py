import sys
import numpy as np

import matplotlib.pyplot as plt

# Get the file name from command line argument
if len(sys.argv) < 2:
    print("Please provide the file name as a command line argument.")
    sys.exit(1)

file_name = sys.argv[1]
figure_file_name = file_name.rsplit(".", 1)[0] + ".png"

# print("figure_file_name:", figure_file_name)
# Read the binary file
data = np.fromfile(file_name, dtype=np.float64)

# Set figure size
plt.figure(figsize=(10, 4))
plt.rcParams.update({'font.size': 16})
# Create scatter plot
plt.scatter(range(len(data)), data, s=1)
# Draw a horizontal line at 1e-10
plt.axhline(y=1e-10, color='r', linestyle='--')
plt.text(0, 1.05e-10, r'$\epsilon=1e^{-10}$', color='r', fontsize=16)
plt.ylabel("Error")
plt.ylim(0,1.5e-10)
plt.savefig(figure_file_name)
plt.show()

print("Max Error:", np.max(data)) 
print("Min Error:", np.min(data))
# Calculate average error
average_error = np.mean(data)
print("Average Error:", average_error)



# Increase font size

