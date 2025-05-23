#!/bin/bash
echo "Start Runing Script"
fixed="../build/examples/fixed/fixed"
adaptive="../build/examples/adaptive/adaptive"
morse_smale="../build/examples/morse_smale/isocontour"
write_vtk="../build/examples/convert/write_vtk"
convert_root_to_vtk="../build/examples/critical_point/convert_root_to_vtk"

derivative_control_point="../build/examples/critical_point/derivative_control_point"
compute_critical_point="../build/examples/critical_point/compute_critical_point"




data_type="sinc_sum"
data_path="../build/examples/${data_type}/"
output_file_prefix="../build/examples/${data_type}/"

mfa_file="${data_path}${data_type}.mfa"
isosurface_file="${data_path}${data_type}_isosurface.dat"
isosurface_file_csv="${data_path}${data_type}_isosurface.csv"
root_file="${data_path}${data_type}_root.dat"
root_file_csv="${data_path}${data_type}_root.csv"


# "${adaptive}" -d 3 -m 2 -q 3 -v 90 -n 200 -i "${data_type}"


# #4.5589 for ackley
# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}.vtk" -m 2 -d 3

# # "${convert_root_to_vtk}" -f "${root_file}" -o "${root_file_csv}" -i "${mfa_file}" -d 0
# # "${convert_root_to_vtk}" -f "${isosurface_file}" -o "${isosurface_file_csv}" -i "${mfa_file}" -d 0

#for schwefel -e 0.48
#for rastrigin -e 0.0012

"${derivative_control_point}" -f "${mfa_file}" -o "${output_file_prefix}${data_type}_cpt.dat"
# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.0012 -t 1e-3
# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.0012 -t 1e-5
# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.0012 -t 1e-9
# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.0012 -t 1e-11
"${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.1 -t 1e-10
"${convert_root_to_vtk}" -f "${output_file_prefix}${data_type}_cp.dat" -o "${output_file_prefix}${data_type}_cp.csv" -i "${mfa_file}" -d 0


