import sys
input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as infile:
    lines = infile.readlines()
    


with open(output_file, 'w') as outfile:
    for line in lines:
        # Check for vertex definition lines (make sure to avoid lines like 'vn' for normals)
        if line.startswith('v '):
            parts = line.strip().split()
            if len(parts) >= 4:
                # parts[0] is "v", parts[1] is x, parts[2] is y, parts[3] is z
                # Set z coordinate to 0 while preserving any additional values (like an optional w)
                new_line = f"{parts[0]} {parts[1]} {parts[2]} 0"
                if len(parts) > 4:
                    new_line += " " + " ".join(parts[4:])
                outfile.write(new_line + "\n")
            else:
                outfile.write(line)
        else:
            # For non-vertex lines, just write them unchanged
            outfile.write(line)




