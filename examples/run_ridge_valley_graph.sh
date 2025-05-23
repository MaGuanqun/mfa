#!/bin/bash
echo "Start Runing Script"
fixed="../build/examples/encode/grid/gridded_2d"
contour="../build/examples/contour/isocontour"
write_vtk="../build/examples/convert/write_vtk"
write_h="../build/examples/convert/write_h"
convert_root_to_vtk="../build/examples/critical_point/convert_root_to_vtk"

derivative_control_point="../build/examples/critical_point/derivative_control_point"
compute_critical_point="../build/examples/critical_point/compute_critical_point"


extract_ttk_js="./python/extract_jacobi_set.py"
set_z_to_zero="./python/set_z_to_zero.py"


extract_block_of_an_obj="./python/extract_block_of_an_obj.py"


count_betti_num="./contour/count_betti_num.py"
show_error_python="./contour/draw_error.py"

# data_type="expotential_function"
# save_folder="expotential"

data_type="cesm"
# data_type="s3d"
# data_type="boussinesq"

save_folder="${data_type}"

extension="ply"

mfa_file="../build/examples/${save_folder}/${data_type}.mfa"
root_file="../build/examples/${save_folder}/${data_type}_root.dat"
root_file_csv="../build/examples/${save_folder}/${data_type}_root.csv"


ridge_valley_file1="../build/examples/${save_folder}/${data_type}_ridge_valley_graph"
edge_type_file="../build/examples/${save_folder}/${data_type}_ridge_valley_graph_edge_type.txt"
ttk_rv_graph_file="../build/examples/${save_folder}/${data_type}_ridge_valley_graph.vtp"
ridge_valley_file_csv1="../build/examples/${save_folder}/${data_type}_ridge_valley_graph.csv"

control_points="../build/examples/${save_folder}/${data_type}_cpt.dat"
critical_point_file="../build/examples/${save_folder}/${data_type}_rt.dat"
critical_point_csv_file="../build/examples/${save_folder}/${data_type}_rt.csv"

error_file="${ridge_valley_file1}_error"

edge_type_file2="../build/examples/${save_folder}/${data_type}_ridge_valley_graph_edge_type_b2.txt"
ttk_rv_graph_file2="../build/examples/${save_folder}/${data_type}_ridge_valley_graph_b2.vtp"

step_size="32"
step_size2="4"
step_size3="8"
step_size4="16"
step_size5="32"
step_size6="64"
step_size7="128"
step_size8="256"

root_finding_epsilon2="1e-6"
root_finding_epsilon3="1e-8"
root_finding_epsilon4="1e-12"

trace_split_grad_square_threshold="1e-4"

connection_threshold="2.0"
connection_threshold2="1.0"
connection_threshold3="1.5"
connection_threshold4="2.5"


if [ "${data_type}" = "cesm" ]; then
    raw_data_file="../build/examples/ori_data/FLDSC_1_1800_3600.dat" 
    root_finding_epsilon="1e-10"
    step_size_half="5e-4"
elif [ "${data_type}" = "s3d" ]; then
    raw_data_file="../build/examples/ori_data/6_small.xyz" 
    root_finding_epsilon="1e-10"
    step_size_half="5e-4"
elif [ "${data_type}" = "expotential_function" ]; then
    root_finding_epsilon="1e-10"
    step_size_half="5e-4"
elif [ "${data_type}" = "boussinesq" ]; then
    raw_data_file="../build/examples/boussinesq/velocity_t2000.bin"  #    
    root_finding_epsilon="1e-10"
    step_size_half="5e-4"
fi

###different blocks for cesm
#for limitations
#675-750-775-850
#2445-2520-1130-1205
#1940-2015-750-825

#cesm scale range: 78-421

# block1="2208-2328-1380-1500" #needs step size = 128



# block1="1473-1593-434-554"

block1="1475.59-1799.51-449.74-773.58"

block1_vtk="82-100-25-43"


# block2="2205-2325-1380-1500"
# block2="1938-2058-1572-1692"

block2="2087.41-2411.34-1367.24-1691.07"
block2_vtk="116-134-76-94"
# block_float1="0.6524-0.6724-0.728-0.7682"

# block1="2452-2562-1171-1281"
# block_float1="0.6524-0.6724-0.728-0.7682"


# block_float2="0.4090-0.4290-0.2535-0.2935"

# block3="2000-2072-897-969"
# block_float3="0.5557-0.5757-0.4986-0.5386"

# b1: 2345-2420-1310-1385 0.651-0.673-0.728-0.770  for write_vtk 0.65157 0.672409 0.728182 0.769872
#b2: 1500-1575-660-735 0.416-0.437-0.366-0.409

# if [ "${data_type}" = "expotential_function" ]; then
# "${fixed}" -d 3 -m 2 -n 301 -v 75 -q 4 -s 0.0 -i "expotential_function" -o "${mfa_file}"
# else 
# "${fixed}" -d 3 -m 2 -f "${raw_data_file}" -i "${data_type}" -q 4 -s 0.0 -o "${mfa_file}"
# fi

# "${derivative_control_point}" -f "${mfa_file}" -o "${control_points}"

# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${control_points}" -o "${critical_point_file}" -e "${step_size_half}" #-s "0.408447-0.429086-0.252918-0.294608"

# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${critical_point_csv_file}" -f "${critical_point_file}" -e "${step_size_half}" 


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size}_b1.obj" -z "${step_size}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block1}" -q "${edge_type_file}"


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size}.obj" -z "${step_size}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k #"${block1}"


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" # -k "${block1}"



# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon2}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon2}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block}""

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon3}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon3}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block}""

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon4}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon4}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block}""

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size3}.obj" -z "${step_size3}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k #"${block1}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size4}.obj" -z "${step_size4}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block1}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size5}.obj" -z "${step_size5}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block1}"





# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size}_b2.obj" -z "${step_size}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block2}" -q "${edge_type_file2}"


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_b2.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block2}"


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size3}_b2.obj" -z "${step_size3}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block2}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size4}_b2.obj" -z "${step_size4}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block2}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size5}_b2.obj" -z "${step_size5}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block2}"



# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size6}.obj" -z "${step_size6}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block2}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size7}.obj" -z "${step_size7}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block2}"

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size8}.obj" -z "${step_size8}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" #-k "${block2}"



# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon2}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon2}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block1}" 

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon3}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon3}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block1}" 

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon4}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon4}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -k "${block1}" 

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${connection_threshold2}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -a "${connection_threshold2}" -k "${block1}" 

# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${connection_threshold3}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -a "${connection_threshold3}"  -k "${block1}" 


# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${connection_threshold4}.obj" -z "${step_size2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}.dat" -a "${connection_threshold4}"  -k "${block1}" 







# "${contour}" -f "${mfa_file}" -r "${root_file}" -b "${ridge_valley_file1}_up.obj" -z "${step_size_2}" -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_file}" -w 0 -s "${error_file}_up.dat" #-k "${block}"





# gdb "${contour}" -ex "set args -f ${mfa_file} -r ${root_file} -b ${ridge_valley_file1}.obj -z ${step_size} -p ${trace_split_grad_square_threshold} -y 0.5 -x ${root_finding_epsilon} -j ${critical_point_file} -w 0 -s ${error_file}.dat"



# "${convert_root_to_vtk}" -f "${ridge_valley_file1}" -o "${ridge_valley_file_csv1}" -i "${mfa_file}" -d 0

# # 

"${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_1.vtk" -m 2 -d 3 -u "${step_size}" -g 1 -z 0

# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_b2.vtk" -m 2 -d 3 -u "${step_size}" -g 1 -z 0 -s "${block2_vtk}" 

# # 
# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_b1.vtk" -m 2 -d 3 -u "${step_size}" -g 1 -z 0 -s "${block1_vtk}" 

# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_up_block1.vtk" -m 2 -d 3 -u "${step_size_2}" -g 1 -z 1 -s "${block_float1}" 

# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_block2.vtk" -m 2 -d 3 -u "${step_size}" -g 1 -z 1 -s "${block_float2}" 

# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_block3.vtk" -m 2 -d 3 -u "${step_size}" -g 1 -z 1 -s "${block_float3}" 

# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_up_block3.vtk" -m 2 -d 3 -u "${step_size_2}" -g 1 -z 1 -s "${block_float3}" 


source ~/enter/etc/profile.d/conda.sh
conda activate mfa
# python ./python/merge_obj_edge_type.py -i "${ridge_valley_file1}${step_size}_b1.obj" -j "${edge_type_file}" -o "${ttk_rv_graph_file}"
# python ./python/merge_obj_edge_type.py -i "${ridge_valley_file1}${step_size}_b2.obj" -j "${edge_type_file2}" -o "${ttk_rv_graph_file2}"

# python "${extract_ttk_js}" "${mfa_file}.vtk" "../build/examples/${save_folder}/rv_ttk_${step_size}_3d.obj" "../build/examples/${save_folder}/rv_ttk_${step_size}.bin"


# python "${extract_ttk_js}" "${mfa_file}_b2.vtk" "../build/examples/${save_folder}/rv_ttk_${step_size}_b2.obj" "../build/examples/${save_folder}/rv_ttk_${step_size}_b2.bin"


# python "${extract_ttk_js}" "${mfa_file}_b1.vtk" "../build/examples/${save_folder}/rv_ttk_${step_size}_b1.obj" "../build/examples/${save_folder}/rv_ttk_${step_size}_b1.bin"

# python "${extract_ttk_js}" "${mfa_file}_up.vtk" "../build/examples/${save_folder}/rv_ttk_up.obj" "../build/examples/${save_folder}/rv_ttk_up.bin"


# python "${extract_ttk_js}" "${mfa_file}_up_block1.vtk" "../build/examples/${save_folder}/rv_ttk_up_block1.obj" "../build/examples/${save_folder}/rv_ttk_up_block1.bin"

# python "${extract_ttk_js}" "${mfa_file}_up_block3.vtk" "../build/examples/${save_folder}/rv_ttk_up_block3.obj" "../build/examples/${save_folder}/rv_ttk_up_block3.bin"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon2}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon4}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size4}.obj"



# python "${count_betti_num}" "${ridge_valley_file1}${step_size5}.obj"






# python "${count_betti_num}" "${ridge_valley_file1}${step_size}_b2.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_b2.obj"


# python "${count_betti_num}" "${ridge_valley_file1}${step_size3}_b2.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size4}_b2.obj"



# python "${count_betti_num}" "${ridge_valley_file1}${step_size5}.obj"


# python "${count_betti_num}" "${ridge_valley_file1}${step_size6}.obj"


# python "${count_betti_num}" "${ridge_valley_file1}${step_size7}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size8}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}_up.obj"
# python "${count_betti_num}" "../build/examples/${save_folder}/rv_ttk_${step_size5}_3d.obj"
# python "${count_betti_num}" "../build/examples/${save_folder}/rv_ttk_${step_size}_3d_b1.obj"


# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon2}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon4}.obj"



# ttk_contour_error_file="../build/examples/contour/compute_function_error"

# "${ttk_contour_error_file}" -f "${mfa_file}" -r "../build/examples/${save_folder}/rv_ttk_${step_size5}.bin" -s "../build/examples/${save_folder}/rv_ttk_${step_size5}.dat" -v 0

# "${ttk_contour_error_file}" -f "${mfa_file}" -r "../build/examples/${save_folder}/rv_ttk_${step_size}_b1.bin" -s "../build/examples/${save_folder}/rv_ttk_${step_size}_b1.dat" -v 0


# python "${set_z_to_zero}" "${ridge_valley_file1}${step_size}_b2.obj" "${ridge_valley_file1}${step_size}_b2_2d.obj"
# python "${set_z_to_zero}" "${ridge_valley_file1}${step_size}_b1.obj" "${ridge_valley_file1}${step_size}_b1_2d.obj"
# python "${set_z_to_zero}" "../build/examples/${save_folder}/rv_ttk_${step_size}.obj" "../build/examples/${save_folder}/rv_ttk_${step_size}_2d.obj"

# python "${set_z_to_zero}" "../build/examples/${save_folder}/rv_ttk_up.obj" "../build/examples/${save_folder}/rv_ttk_up_2d.obj"

# python "${extract_block_of_an_obj}" "${ridge_valley_file1}_2d.obj" "${ridge_valley_file1}_2d_block3.obj" "${block3}"
# python "${extract_block_of_an_obj}" "${ridge_valley_file1}_up_2d.obj" "${ridge_valley_file1}_up_2d_block3.obj" "${block3}"

# python "${extract_block_of_an_obj}" "../build/examples/${save_folder}/rv_ttk_${step_size}_2d.obj" "../build/examples/${save_folder}/rv_ttk_${step_size}_2d_block3.obj" "${block3}"

# python "${show_error_python}" "../build/examples/${save_folder}/rv_ttk_${step_size}.dat"