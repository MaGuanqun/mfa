import pandas as pd
import argparse

def extract_points(range,path,input_csv,data,output_csv):
# Step 1: Read the CSV file
    df = pd.read_csv(path+input_csv+'.csv')
    
    if data[:3]=='rti':
        data_size=[143,255,255]
    elif data[:3]=='qmc':
        data_size=[68,68,114]
    
    # Step 2: Define the cube's bounds
    # Let's say the cube is defined from (min_x, min_y, min_z) to (max_x, max_y, max_z)
    min_x=data_size[0]*range[0]
    min_y=data_size[1]*range[2]
    min_z=data_size[2]*range[4]
    max_x=data_size[0]*range[1]
    max_y=data_size[1]*range[3]
    max_z=data_size[2]*range[5]
    
    if input_csv[:3]=='ttk':            
        # Step 3: Filter points within the cube
        filtered_df = df[(df['PositionX'] >= min_x) & (df['PositionX'] <= max_x) &
                        (df['PositionY'] >= min_y) & (df['PositionY'] <= max_y) &
                        (df['PositionZ'] >= min_z) & (df['PositionZ'] <= max_z)]
    else:
        filtered_df = df[(df['x0'] >= min_x) & (df['x0'] <= max_x) &
                        (df['x1'] >= min_y) & (df['x1'] <= max_y) &
                        (df['x2'] >= min_z) & (df['x2'] <= max_z)]
    # Step 4: Save the filtered points to a new CSV file
    print(path+output_csv+".csv")
    filtered_df.to_csv(path+output_csv+".csv", index=False)



parser = argparse.ArgumentParser(description='plot critical points.')

parser.add_argument('-p', '--path', type=str, default='path', help='path to all files')
parser.add_argument('-f', '--input_csv', type=str, default='cesm_1', help='cesm_1, cesm_2...')
parser.add_argument('-o', '--output_csv', type=str, default='ttk_1', help='ttk critical point file name')
parser.add_argument('-r', '--local_range', type=str, default='0,1,0,1,0,1', help='the range of points')
parser.add_argument('-d', '--data_name', type=str, default='cesm', help='cesm, qmcpack..')

args = parser.parse_args()

range =  tuple(map(float, args.local_range.split('-')))
path=args.path
input_csv=args.input_csv
output_csv=args.output_csv
data=args.data_name

extract_points(range,path,input_csv,data,output_csv)

