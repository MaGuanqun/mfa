import csv

def extract_points_in_block(input_file, output_file, block):
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip the header row
        points = []
        for row in reader:
            x, y = float(row[0]), float(row[1])
            if block[0] <= x <= block[2] and block[1] <= y <= block[3]:
                points.append(row)

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(points)

# Example usage
input_file = '../../build/examples/cesm/cesm_rt.csv'
output_file1 = '../../build/examples/cesm/cesm_rt_b1.csv'
output_file2= '../../build/examples/cesm/cesm_rt_b2.csv'

block1 = [1475.5900000000001, 449.74999999999994,1799.5,  773.5699999999999]  # Specify the block as [x_min, y_min, x_max, y_max]
block2=[2087.42,1367.2399999999998, 2411.33,  1691.06]

extract_points_in_block(input_file, output_file1, block1)
extract_points_in_block(input_file, output_file2, block2)