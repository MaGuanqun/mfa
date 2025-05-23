import pandas as pd
import numpy as np
from scipy.spatial import KDTree
import argparse

def read_csv_to_3d_array(file_path):
    """Read CSV file and convert to a numpy array of 3D points."""
    df = pd.read_csv(file_path)  # Adjust if your CSV has a header
    return df[['x0', 'x1', 'x2', 'x3', 'x4']].to_numpy()


def read_ttk_csv_to_3d_array(file_path):
    """Read CSV file and convert to a numpy array of 3D points."""
    df = pd.read_csv(file_path)  # Adjust if your CSV has a header
    return df[['PositionX', 'PositionY', 'PositionZ', 'CriticalType']].to_numpy()


def read_csv_to_2d_array(file_path):
    """Read CSV file and convert to a numpy array of 3D points."""
    df = pd.read_csv(file_path)  # Adjust if your CSV has a header
    return df[['x0', 'x1', 'x2', 'x3']].to_numpy()

def read_ttk_csv_to_2d_array(file_path):
    """Read CSV file and convert to a numpy array of 3D points."""
    df = pd.read_csv(file_path)  # Adjust if your CSV has a header
    return df[['PositionX', 'PositionY','PositionZ', 'CriticalType']].to_numpy()

parser = argparse.ArgumentParser(description='find overlap.')

parser.add_argument('-p', '--path', type=str, default='path', help='path to all files')
parser.add_argument('-d', '--threshold', type=str, default='1.0', help='threshold...')
parser.add_argument('-t','--ttk_csv', type=str, default='default', help='ttk critical point file name')
parser.add_argument('-o','--output', type=str, default='default', help='output csv name')
parser.add_argument('-f', '--cpe_csv', type=str, default='default', help='cpe-mfa critical point file name')
parser.add_argument('-g', '--extension', type=str, default='ply', help='the extrension of original file to determine dimension')
parser.add_argument('-o1c','--output_array1_close', type=str, default='default', help='output csv name of points in array1 that have close point in array2')
parser.add_argument('-o2c','--output_array2_close', type=str, default='default', help='output csv name of points in array2 that have close point in array1')
parser.add_argument('-o1f','--output_array1_far', type=str, default='default', help='output csv name of points in array1 that dont have close point in array2')
parser.add_argument('-o2f','--output_array2_far', type=str, default='default', help='output csv name of points in array2 that dont have close point in array1')

args = parser.parse_args()


k = float(args.threshold)  # Example distance, change according to your requirement



if args.extension=="ply":
    dimension=2
else:
    dimension=3


input_cpe_file_name=args.path+args.cpe_csv+".csv"
input_ttk_file_name=args.path+args.ttk_csv+".csv"
output_file_name=args.path+args.output+".csv"

output_array1_c_file_name = args.path+args.output_array1_close+".csv"
output_array2_c_file_name = args.path+args.output_array2_close+".csv"
output_array1_f_file_name = args.path+args.output_array1_far+".csv"
output_array2_f_file_name = args.path+args.output_array2_far+".csv"


# Reading and converting CSV files to 3D point arrays

if dimension == 3:
    array1_all = read_csv_to_3d_array(input_cpe_file_name)
    array2_all = read_ttk_csv_to_3d_array(input_ttk_file_name)
    array1=array1_all[:,:3]
    array2=array2_all[:,:3]
else: 
    array1_all = read_csv_to_2d_array(input_cpe_file_name)
    array2_all = read_ttk_csv_to_2d_array(input_ttk_file_name)
    array1=array1_all[:,:2]
    array2=array2_all[:,:2]

# Construct KD-tree for the second array of points
tree = KDTree(array2)
points_with_neighbors =[]
labeled_points_array1 = set()
labeled_points_array2  = set()

labeled_points_index_same = []

for index, point in enumerate(array1):
    neighbors_index = tree.query_ball_point(point, k)
    if neighbors_index:
        labeled_points_array1.add(index)
        points_with_neighbors.append(point)
        for neighbor_index in neighbors_index:
            labeled_points_array2.add(neighbor_index)
        temp=int(array1_all[index][-1])
        if dimension==2:
            if temp == 2:
                temp = 3
        for neighbor_index in neighbors_index:          
            if temp==int(array2_all[neighbor_index][-1]):
                labeled_points_index_same.append(point)
                break
            
            

# Convert the list of points to a DataFrame

# if dimension==3 :
#     points_with_neighbors_df = pd.DataFrame(points_with_neighbors, columns=['x0', 'x1', 'x2'])
# else:
#     points_with_neighbors_df = pd.DataFrame(points_with_neighbors, columns=['x0', 'x1'])

# points_with_neighbors_df.to_csv(output_file_name, index=False)

print("distance "+ str(k)+ " CPE "+str(len(array1))+" TTK "+str(len(array2))+ " Overlap "+str(len(points_with_neighbors)))
print("index matching "+str(len(labeled_points_index_same)))


labeled_array2 = [array2_all[i] for i in labeled_points_array2]
labeled_array1 = [array1_all[i] for i in labeled_points_array1]

unlabeled_array2 = [array2_all[i] for i in range(len(array2)) if i not in labeled_points_array2]
unlabeled_array1 = [array1_all[i] for i in range(len(array1)) if i not in labeled_points_array1]

if dimension==3 :
    points_with_neighbors_array1 = pd.DataFrame(labeled_array1, columns=['x0', 'x1', 'x2', 'x3','x4'])
    points_with_neighbors_array2 = pd.DataFrame(labeled_array2, columns=['x0', 'x1', 'x2', 'x3'])
    points_without_neighbors_array1 = pd.DataFrame(unlabeled_array1, columns=['x0', 'x1', 'x2', 'x3','x4'])
    points_without_neighbors_array2 = pd.DataFrame(unlabeled_array2, columns=['x0', 'x1', 'x2', 'x3'])
else:
    points_with_neighbors_array1 = pd.DataFrame(labeled_array1, columns=['x0', 'x1', 'x2','x3'])
    points_with_neighbors_array2 = pd.DataFrame(labeled_array2, columns=['x0', 'x1', 'x2','x3'])
    points_without_neighbors_array1 = pd.DataFrame(unlabeled_array1, columns=['x0', 'x1', 'x2','x3'])
    points_without_neighbors_array2 = pd.DataFrame(unlabeled_array2, columns=['x0', 'x1', 'x2','x3'])

points_with_neighbors_array1.to_csv(output_array1_c_file_name, index=False)
points_with_neighbors_array2.to_csv(output_array2_c_file_name, index=False)
points_without_neighbors_array1.to_csv(output_array1_f_file_name, index=False)
points_without_neighbors_array2.to_csv(output_array2_f_file_name, index=False)

# Print the number of points that have a neighbor within distance k
# print(f"Number of overlap closer than {k}: {len(points_with_neighbors)}")

