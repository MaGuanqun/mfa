#!/bin/bash
echo "Start Runing Script"
fixed="../build/examples/encode/grid/gridded_2d"
contour="../build/examples/contour/isocontour"
write_vtk="../build/examples/convert/write_vtk"
convert_root_to_vtk="../build/examples/critical_point/convert_root_to_vtk"


derivative_control_point="../build/examples/critical_point/derivative_control_point"
compute_critical_point="../build/examples/critical_point/compute_critical_point"

show_error_python="./contour/draw_error.py"
count_betti_num="./contour/count_betti_num.py"

extract_ttk_contour="./python/extract_contour.py"

set_z_to_zero="./python/set_z_to_zero.py"

# data_type="schwefel" #"sinc_sum"
data_type="sinc_sum"
# data_type="expotential_function"
# save_folder="expotential"

sinc_origin_data="../build/examples/ori_func_sample/sinc_sum_2d_triangles.vtk"

schwefel_origin_data="../build/examples/ori_func_sample/schwefel_2d_triangles.vtk"

# data_type="cesm"
# data_type="s3d"
save_folder="${data_type}"

extension="ply"

mfa_file="../build/examples/${save_folder}/${data_type}.mfa"
root_file="../build/examples/${save_folder}/${data_type}_root.dat"
root_file_csv="../build/examples/${save_folder}/${data_type}_root.csv"


ridge_valley_file1="../build/examples/${save_folder}/${data_type}_isocontour"
ridge_valley_file_csv1="../build/examples/${save_folder}/${data_type}_isocontour.csv"

error_file="${ridge_valley_file1}_error"


control_points="../build/examples/${save_folder}/${data_type}_cpt.dat"
critical_point_file="../build/examples/${save_folder}/${data_type}_rt.dat"
critical_point_csv_file="../build/examples/${save_folder}/${data_type}_rt.csv"


trace_split_grad_square_threshold="1e-4"
step_size="2"
step_size2="4"
step_size3="8"
step_size4="16"
step_size5="32"
step_size6="64"
step_size7="128"
step_size8="256"

root_finding_epsilon="1e-10"

root_finding_epsilon2="1e-6"
root_finding_epsilon3="1e-8"
root_finding_epsilon4="1e-12"

connection_threshold="2.0"

connection_threshold2="1.0"
connection_threshold3="1.5"
connection_threshold4="2.5"



if [ "${data_type}" = "cesm" ]; then
    raw_data_file="../build/examples/ori_data/FLDSC_1_1800_3600.dat" 
elif [ "${data_type}" = "s3d" ]; then
    value1="30" #50,60
    value2="60"
    value3="60"
    raw_data_file="../build/examples/ori_data/6_small.xyz" 
elif [ "${data_type}" = "schwefel" ]; then
    value1="100" #500
    value2="100"
elif [ "${data_type}" = "sinc_sum" ]; then
    value1="0.33" #0.73
    value2="0.79"
elif [ "${data_type}" = "expotential_function" ]; then
    value1="1.05"
fi



# if [ "${data_type}" = "schwefel" ]; then
# "${fixed}" -d 3 -m 2 -n 301 -v 75 -q 4 -s 0.0 -i "schwefel" -o "${mfa_file}"
# else 
# "${fixed}" -d 3 -m 2 -f "${raw_data_file}" -i "${data_type}" -q 4 -s 0.0 -o "${mfa_file}"
# fi

# "${derivative_control_point}" -f "${mfa_file}" -o "${control_points}"

# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${control_points}" -o "${critical_point_file}" -e "${step_size}"

# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${critical_point_csv_file}" -f "${critical_point_file}" -e "${step_size}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value1}_${step_size}.obj" -z "${step_size}" -v "${value1}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value1}.dat"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size}.obj" -z "${step_size}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2


"${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}_${connection_threshold2}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a "${connection_threshold2}"

"${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}_${connection_threshold3}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a "${connection_threshold3}"

"${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}_${connection_threshold4}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a "${connection_threshold4}"



# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}_${root_finding_epsilon2}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon2}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}_${root_finding_epsilon3}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon3}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size2}_${root_finding_epsilon4}.obj" -z "${step_size2}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon4}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" 

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size3}.obj" -z "${step_size3}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size4}.obj" -z "${step_size4}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size5}.obj" -z "${step_size5}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size6}.obj" -z "${step_size6}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size7}.obj" -z "${step_size7}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value2}_${step_size8}.obj" -z "${step_size8}" -v "${value2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 0 -s "${error_file}_${value2}.dat" -a 2


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_${value3}_${step_size}.obj" -z "${step_size}" -v "${value3}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -u 1 -w 1 -s "${error_file}_${value3}.dat"







# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}.vtk" -m 2 -d 3 -u "${step_size}" -g 0 -z 0 
# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_${step_size}.vtk" -m 2 -d 3 -u "${step_size}" -g 0 -z 1


source ~/enter/etc/profile.d/conda.sh
conda activate mfa

# # python "${show_error_python}" "${error_file}_${value1}.dat"


# python "${extract_ttk_contour}" "${mfa_file}.vtk" "${value1}" "../build/examples/${save_folder}/contour_${value1}_ttk_3d.obj" "../build/examples/${save_folder}/contour_${value1}_ttk.bin"

# python "${extract_ttk_contour}" "${mfa_file}.vtk" "${value2}" "../build/examples/${save_folder}/contour_${value2}_ttk_3d.obj" "../build/examples/${save_folder}/contour_${value2}_ttk.bin"

# python "${extract_ttk_contour}" "${mfa_file}.vtk" "${value3}" "../build/examples/${save_folder}/contour_${value3}_ttk_3d.obj" "../build/examples/${save_folder}/contour_${value3}_ttk.bin"



# python "${extract_ttk_contour}" "${mfa_file}_${step_size}.vtk" "${value1}" "../build/examples/${save_folder}/contour_${value1}_ttk.obj" "../build/examples/${save_folder}/contour_${value1}_ttk.bin"

# python "${extract_ttk_contour}" "${mfa_file}_${step_size}.vtk" "${value2}" "../build/examples/${save_folder}/contour_${value2}_ttk.obj" "../build/examples/${save_folder}/contour_${value2}_ttk.bin"

# python "${extract_ttk_contour}" "${mfa_file}_${step_size}.vtk" "${value3}" "../build/examples/${save_folder}/contour_${value3}_ttk.obj" "../build/examples/${save_folder}/contour_${value3}_ttk.bin"




# python "${extract_ttk_contour}" "${mfa_file}_3d_up.vtk" "${value1}" "../build/examples/${save_folder}/contour_${value1}_ttk_up.obj" "../build/examples/${save_folder}/contour_${value1}_ttk_up.bin"
# python "${extract_ttk_contour}" "${mfa_file}_3d_${step_size}.vtk" "${value2}" "../build/examples/${save_folder}/contour_${value2}_ttk_${step_size}.obj" "../build/examples/${save_folder}/contour_${value2}_ttk_${step_size}.bin"

# python "${extract_ttk_contour}" "${mfa_file}_3d_up.vtk" "${value3}" "../build/examples/${save_folder}/contour_${value3}_ttk_up.obj" "../build/examples/${save_folder}/contour_${value3}_ttk_up.bin"

# python "${extract_ttk_contour}" "${sinc_origin_data}" "${value1}" "../build/examples/ori_func_sample/sinc_ori_value1.obj" "../build/examples/ori_func_sample/sinc_ori_value1.bin"
# python "${extract_ttk_contour}" "${sinc_origin_data}" "${value2}" "../build/examples/ori_func_sample/sinc_ori_value2.obj" "../build/examples/ori_func_sample/sinc_ori_value2.bin"


# python "${count_betti_num}" "${ridge_valley_file1}_${value1}_${step_size}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_value.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value3}_value.obj"

# python "${count_betti_num}" "${ridge_valley_file1}_30_up.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}.obj"
python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}_${connection_threshold2}.obj"
python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}_${connection_threshold3}.obj"
python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}_${connection_threshold4}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}_${root_finding_epsilon2}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}_${root_finding_epsilon3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size2}_${root_finding_epsilon4}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size4}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size5}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size6}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size7}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_${value2}_${step_size8}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}_60_up.obj"



# python "${count_betti_num}" "../build/examples/ori_func_sample/sinc_ori_value1.obj"
# python "${count_betti_num}" "../build/examples/ori_func_sample/sinc_ori_value2.obj"
# python "${count_betti_num}" "../build/examples/${save_folder}/contour_${value1}_ttk.obj" 
# python "${count_betti_num}" "../build/examples/${save_folder}/contour_${value2}_ttk.obj" 
# python "${count_betti_num}" "../build/examples/${save_folder}/contour_${value3}_ttk.obj" 


# python "${set_z_to_zero}" "${ridge_valley_file1}_${value1}_${step_size}.obj" "${ridge_valley_file1}_${value1}_${step_size}_2d.obj"
# python "${set_z_to_zero}" "${ridge_valley_file1}_${value2}_${step_size}.obj" "${ridge_valley_file1}_${value2}_${step_size}_2d.obj"
# python "${set_z_to_zero}" "${ridge_valley_file1}_${value3}_${step_size}.obj" "${ridge_valley_file1}_${value3}_${step_size}_2d.obj"




# "${convert_root_to_vtk}" -f "${root_file}" -o "${root_file_csv}" -i "${mfa_file}" -d 0 -u 5


ttk_contour_error_file="../build/examples/contour/compute_function_error"

# "${ttk_contour_error_file}" -f "${mfa_file}" -r "../build/examples/${save_folder}/contour_${value1}_ttk.bin" -s "${error_file}_${value1}_contour.dat" -v "${value1}"
# "${ttk_contour_error_file}" -f "${mfa_file}" -r "../build/examples/${save_folder}/contour_${value2}_ttk.bin" -s "${error_file}_${value2}_contour.dat" -v "${value2}"
# "${ttk_contour_error_file}" -f "${mfa_file}" -r "../build/examples/${save_folder}/contour_${value3}_ttk.bin" -s "${error_file}_${value3}_contour.dat" -v "${value3}"