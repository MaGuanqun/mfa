#!/bin/bash
echo "Start Runing Script"
analytical="../build/examples/encode/analytical/analytical"
gridded_2d="../build/examples/encode/grid/gridded_2d"

derivative_control_point="../build/examples/critical_point/derivative_control_point"
compute_critical_point="../build/examples/critical_point/compute_critical_point"
convert_root_to_vtk="../build/examples/critical_point/convert_root_to_vtk"
convert_cpt_to_vtk="../build/examples/critical_point/convert_cpt_to_vtk"


data_type="cesm"
data_path="../build/examples/critical_point/"
output_file_prefix="../build/examples/${data_type}/"

mfa_file="${data_path}${data_type}.mfa"

raw_data_file="../build/examples/ori_data/FLDSC_1_1800_3600.dat" #rti #dd07g_xxsmall_le

# "${analytical}" -d 3 -m 2 -q 3 -v 90 -n 200 -i "${data_type}" -f "${mfa_file}"

"${gridded_2d}" -d 3 -q 2 -f "${raw_data_file}" -i "${data_type}" -o "${mfa_file}" -a 1 -e "1.0e-2"



"${derivative_control_point}" -f "${mfa_file}" -o "${output_file_prefix}${data_type}_cpt.dat"
# "${convert_cpt_to_vtk}" -f "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_gradient_control"
"${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.1 -t 1e-10
# "${convert_root_to_vtk}" -f "${output_file_prefix}${data_type}_cp.dat" -o "${output_file_prefix}${data_type}_cp.csv" -i "${mfa_file}" -d 0


