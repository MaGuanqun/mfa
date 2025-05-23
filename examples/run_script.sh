#!/bin/bash
echo "Start Runing Script"
convert_root_to_vtk="../build/examples/critical_point/convert_root_to_vtk"
write_vtk="../build/examples/convert/write_vtk"
derivative_control_point="../build/examples/critical_point/derivative_control_point"
compute_critical_point="../build/examples/critical_point/compute_critical_point"

data_type="rti"
extension="vtk"  #ply or vtk
if [ "${data_type}" = "cesm" ]; then
    anchor_pos="1.0,1.0" #cesm
    figure_size="12,6" #cesm  
    block_1="0-0.2-0.55-0.75"
    block_2="0.25-0.45-0.25-0.45"
    block_3="0.6-0.8-0.6-0.8"
    block_4="0.05-0.06-0.704-0.714"
    block_5="0.077-0.081-0.723-0.733"
    block_10="0.2-0.4-0.55-0.75"
    block_11="0.75-0.95-0.48-0.6804"
    adaptive="../build/examples/encode/grid/gridded_2d"
    # figure_size="6,7" #cesm
elif [ "${data_type}" = "rti" ]; then
    anchor_pos="0.6,0.8" #rti
    figure_size="6,6" #rti
    block_1="0.52-0.6293-0.259-0.307-0.345-0.385"
    block_2="0.4-0.51-0.09-0.1373-0.61-0.65"
    block_3="0.31-0.42-0.2-0.247-0.80-0.84"
    raw_data_file="../build/examples/ori_data/rti_144_256_256.xyz" #rti #dd07g_xxsmall_le
    block_1_local="0.526-0.5535-0.29337-0.30587-0.368-0.378"
    block_2_local="0.4-0.4275-0.0902-0.1018-0.6118-0.6218"
    block_3_local="0.31-0.3375-0.2-0.2125-0.81-0.82"
    block_4="0.522-0.532-0.297-0.307-0.369-0.379" # have to round to 0.1 instead of 1 when setting the range
    block_5="0.399-0.409-0.110-0.120-0.642-0.652" #block 2
    block_6="0.319-0.329-0.238-0.248-0.81-0.82" #block 3
    block_10="0.4-0.51-0.09-0.1373-0.20-0.24"
    block_11="0.35-0.46-0.15-0.197-0.3-0.34"
    adaptive="../build/examples/encode/grid/gridded_3d"
    # anchor_pos="0.6,0.9" #rti
    

    # 6&8 in block 2, 7 in block 3 
    #  block_4="0.5944055944055944-0.6083916083916084-0.25960784313725493-0.26745098039215687-0.3501960784313725-0.3580392156862745" #block 1
    # block_7="0.4125874125874126-0.42657342657342656-0.23921568627450981-0.24705882352941178-0.8-0.807843137254902" #block2
    # 
elif [ "${data_type}" = "nek" ]; then
    anchor_pos="1.2,1.0" #nek
    figure_size="8,6" #nek
    block_1="0.6-0.9-0.35-0.653"
    block_2="0.0-0.3-0.0-0.3"
    block_3="0.35-0.653-0.7-1.0"
    raw_data_file="../build/examples/ori_data/0.xyz"
    adaptive="../build/examples/encode/grid/gridded_2d"
elif [ "${data_type}" = "qmcpack" ]; then
    echo "select qmcpack"
    anchor_pos="0.8,0.8" #qmcpack
    figure_size="10,10" #qmcpack
    block_1="0.7-1.0-0.7-1.0-0.4-0.7"
    block_2="0.35-0.65-0.35-0.65-0.6-0.895"
    block_3="0.6-0.9-0.1-0.4-0.3-0.6"
    block_1_local="0.82-0.93-0.82-0.93-0.42-0.53"
    block_2_local="0.44-0.55-0.44-0.55-0.69-0.8"
    block_3_local="0.75-0.86-0.1-0.21-0.42-0.53"
    block_4="0.763235294117647-0.8073529411764706-0.7058823529411765-0.75-0.4473684210526316-0.4614035087719298" #block 1
    block_5="0.43823529411764706-0.4838235294117647-0.43676470588235294-0.48088235294117654-0.7780701754385965-0.7921052631578948"
    block_6="0.8529411764705882-0.8970588235294118-0.16911764705882354-0.21323529411764705-0.5096491228070176-0.5236842105263159"
    block_10="0.35-0.65-0.35-0.65-0.3-0.6"
    block_11="0.65-0.941-0.4-0.691-0.35-0.65"
    figure_size="7,8" #qmcpack
    anchor_pos="0.3,1.05" #qmcpack
    adaptive="../build/examples/encode/grid/gridded_3d"
    raw_data_file="../build/examples/ori_data/qmcpack_288_115_69_69.pre.f32" #rti
elif [ "${data_type}" = "s3d" ]; then
    anchor_pos="1.2,1.0" #s3d
    figure_size="10,6" #s3d
    block_1="0.6-0.9-0.5-0.8"
    block_2="0.25-0.55-0.35-0.65"
    block_3="0.1-0.4-0.052-0.35"
    block_5="0.28-0.31-0.365-0.395"
    block_6="0.225-0.255-0.185-0.215"
    block_7="0.225-0.236-0.108-0.119"
    block_8="0.307-0.317-0.099-0.106"
    block_9="0.281-0.288-0.377-0.388"
    block_10="0.4-0.7-0.052-0.35"
    block_11="0.7-1.0-0.25-0.55"
    raw_data_file="../build/examples/ori_data/6_small.xyz"
    adaptive="../build/examples/encode/grid/gridded_2d"
    # figure_size="6,6" #s3d
    # anchor_pos="1.25,1.0" #s3d
    # figure_size="6,7" #s3d
else
    anchor_pos="1.2,1.0" #default
    figure_size="8,6" #default
fi
path="../build/examples/${data_type}/"
output_file_prefix="${path}${data_type}_"
mfa_file="../build/examples/${data_type}/${data_type}.mfa"
mfa_file_test="../build/examples/${data_type}/${data_type}_test.mfa"
root_file="${output_file_prefix}rt.dat"
root_file_up="${output_file_prefix}rt_0_1.dat"
root_file_up_2="${output_file_prefix}rt_0_0_1.dat"
root_temp_Save="${output_file_prefix}rt_temp.dat"


# error_2="0.1998"
error_2="0.0999"
up_ratio="10"
error_1="0.999"
up_ratio_2="100"
error_3="0.009999"
# error_3="0.0024975"


# #approx mfa
# "${adaptive}" -e "1.0e-2" -d "4" -q "2" -i "$data_type" -a "1" -f "${raw_data_file}" -o "${mfa_file_test}"   

# # # # # #compute critical points
# "${derivative_control_point}" -f "${mfa_file}" -o "${output_file_prefix}cpt.dat"
# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file}" -e "${error_1}"
"${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up}" -e "${error_2}"
# "${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up_2}" -e "${error_3}"


# "${compute_critical_point}" -s "${block_1}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file}" -e "${error_1}"
# "${compute_critical_point}" -s "${block_2}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file}" -e "${error_1}"
# "${compute_critical_point}" -s "${block_3}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file}" -e "${error_1}"

# "${compute_critical_point}" -s "${block_1}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up}" -e "${error_2}"
# "${compute_critical_point}" -s "${block_2}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up}" -e "${error_2}"
# "${compute_critical_point}" -s "${block_3}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up}" -e "${error_2}"

# "${compute_critical_point}" -s "${block_11}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up}" -e "${error_2}"

# "${compute_critical_point}" -s "${block_4}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up_2}" -e "${error_3}"
# "${compute_critical_point}" -s "${block_5}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up_2}" -e "${error_3}"
# "${compute_critical_point}" -s "${block_6}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}cpt.dat" -o "${root_file_up_2}" -e "${error_3}"

# # # # # # # # extract root from .dat file
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${path}${data_type}.csv" -f "${root_file}" -e "${error_1}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}1.csv" -s "${block_1}" -f "${root_file}" -e "${error_1}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}2.csv" -s "${block_2}" -f "${root_file}" -e "${error_1}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}3.csv" -s "${block_3}" -f "${root_file}" -e "${error_1}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${path}${data_type}_up.csv" -f "${root_file_up}" -e "${error_2}"
"${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}1_up.csv" -s "${block_1}" -f "${root_file_up}" -e "${error_2}"
"${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}2_up.csv" -s "${block_2}" -f "${root_file_up}" -e "${error_2}"
"${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}3_up.csv" -s "${block_3}" -f "${root_file_up}" -e "${error_2}"

# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}11_up.csv" -s "${block_11}" -f "${root_file_up}" -e "${error_2}"

# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}4_up.csv" -s "${block_4}" -f "${root_file_up_2}" -e "${error_3}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}5_up.csv" -s "${block_5}" -f "${root_file_up_2}" -e "${error_3}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}6_up.csv" -s "${block_6}" -f "${root_file_up_2}" -e "${error_3}"

# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}8_up.csv" -s "${block_8}" -f "${root_file_up_2}" -e "${error_3}"
#  "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}9_up.csv" -s "${block_9}" -f "${root_file_up_2}" -e "${error_3}"

# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}10.csv" -s "${block_10}" -f "${root_file}" -e "${error_1}"
"${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}10_up.csv" -s "${block_10}" -f "${root_file_up}" -e "${error_2}"
# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}11.csv" -s "${block_11}" -f "${root_file}" -e "${error_1}"
"${convert_root_to_vtk}" -i "${mfa_file}" -o "${output_file_prefix}11_up.csv" -s "${block_11}" -f "${root_file_up}" -e "${error_2}"


# if [ "${extension}" = "vtk" ]; then
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -t "${path}${data_type}"
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_1}" -t "${output_file_prefix}1"
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_2}" -t "${output_file_prefix}2"
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_3}" -t "${output_file_prefix}3"
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_1}" -t "${output_file_prefix}1_up"
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_2}" -t "${output_file_prefix}2_up"
#     "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_3}" -t "${output_file_prefix}3_up"

# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_4}" -t "${output_file_prefix}4"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_5}" -t "${output_file_prefix}5"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_6}" -t "${output_file_prefix}6"

# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_4}" -t "${output_file_prefix}4_up"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_5}" -t "${output_file_prefix}5_up"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_6}" -t "${output_file_prefix}6_up"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_7}" -t "${output_file_prefix}7_up"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_8}" -t "${output_file_prefix}8_up"

# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_10}" -t "${output_file_prefix}10"
# # "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "10" -s "${block_10}" -t "${output_file_prefix}10_up"

# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_11}" -t "${output_file_prefix}11"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "10" -s "${block_11}" -t "${output_file_prefix}11_up"

# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_12}" -t "${output_file_prefix}12"
# "${write_vtk}" -m "3" -d "4" -f "${mfa_file}" -i "${data_type}" -u "10" -s "${block_12}" -t "${output_file_prefix}12_up"

#  else # 2d domain
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "1" -n "${path}${data_type}.ply" -s "0-1-0-1"
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_1}" -n "${output_file_prefix}1.ply"
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_2}" -n "${output_file_prefix}2.ply"
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_3}" -n "${output_file_prefix}3.ply"
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_1}" -n "${output_file_prefix}1_up.ply"
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_2}" -n "${output_file_prefix}2_up.ply"
    # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio}" -s "${block_3}" -n "${output_file_prefix}3_up.ply"
# # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_4}" -n "${output_file_prefix}4_up.ply"
# "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_5}" -n "${output_file_prefix}5_up.ply"
# "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_6}" -n "${output_file_prefix}6_up.ply"
# # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_7}" -n "${output_file_prefix}7_up.ply"
# # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_8}" -n "${output_file_prefix}8_up.ply"
# # "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "${up_ratio_2}" -s "${block_9}" -n "${output_file_prefix}9_up.ply"

# "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_10}" -n "${output_file_prefix}10.ply"
# "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "10" -s "${block_10}" -n "${output_file_prefix}10_up.ply"

# "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "1" -s "${block_11}" -n "${output_file_prefix}11.ply"
# "${write_vtk}" -m "2" -d "3" -f "${mfa_file}" -i "${data_type}" -u "10" -s "${block_11}" -n "${output_file_prefix}11_up.ply"


# # fi


# python to create critical point
source ~/enter/etc/profile.d/conda.sh
conda activate mfa

# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}.${extension}" -o "ttk.csv"

# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_1.${extension}" -o "ttk_1.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_2.${extension}" -o "ttk_2.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_3.${extension}" -o "ttk_3.csv"

# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_1_up.${extension}" -o "ttk_1_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_2_up.${extension}" -o "ttk_2_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_3_up.${extension}" -o "ttk_3_up.csv"

# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_4_up.${extension}" -o "ttk_4_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_5_up.${extension}" -o "ttk_5_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_6_up.${extension}" -o "ttk_6_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_7_up.${extension}" -o "ttk_7_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_8_up.${extension}" -o "ttk_8_up.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_9_up.${extension}" -o "ttk_9_up.csv"

# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_10.${extension}" -o "ttk_10.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_10_up.${extension}" -o "ttk_10_up.csv"
# 
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_11.${extension}" -o "ttk_11.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_11_up.${extension}" -o "ttk_11_up.csv"

# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_12.${extension}" -o "ttk_12.csv"
# python ./python/ttk_critical_point.py -p "${path}" -i "${data_type}_12_up.${extension}" -o "ttk_12_up.csv"

source ~/enter/etc/profile.d/conda.sh
conda activate ttk


# python ./python/extract_local.py -p "${path}" -f "${data_type}_1" -o "${data_type}-1-local" -r "${block_1_local}" -d "${data_type}"
# python ./python/extract_local.py -p "${path}" -f "ttk_1" -o "ttk-1-local" -r "${block_1_local}" -d "${data_type}"
# python ./python/extract_local.py -p "${path}" -f "${data_type}_1_up" -o "${data_type}-1-up-local" -r "${block_1_local}" -d "${data_type}"
# python ./python/extract_local.py -p "${path}" -f "ttk_1_up" -o "ttk-1-up-local" -r "${block_1_local}" -d "${data_type}"

# python ./python/extract_local.py -p "${path}" -f "${data_type}_2" -o "${data_type}-2-local" -r "${block_2_local}" -d "${data_type}"
# python ./python/extract_local.py -p "${path}" -f "ttk_2" -o "ttk-2-local" -r "${block_2_local}" -d "${data_type}"
# python ./python/extract_local.py -p "${path}" -f "${data_type}_2_up" -o "${data_type}-2-up-local" -r "${block_2_local}" -d "${data_type}"
# python ./python/extract_local.py -p "${path}" -f "ttk_2_up" -o "ttk-2-up-local" -r "${block_2_local}" -d "${data_type}"

# python ./python/extract_local.py -p "${path}" -f "${data_type}_3" -o "${data_type}-3-local" -r "${block_3_local}" -d "${data_type}"

# python ./python/extract_local.py -p "${path}" -f "ttk_3" -o "ttk-3-local" -r "${block_3_local}" -d "${data_type}"

# python ./python/extract_local.py -p "${path}" -f "${data_type}_3_up" -o "${data_type}-3-up-local" -r "${block_3_local}" -d "${data_type}"

# python ./python/extract_local.py -p "${path}" -f "ttk_3_up" -o "ttk-3-up-local" -r "${block_3_local}" -d "${data_type}"



# if [ "${extension}" = "vtk" ]; then
#     echo "The extension is vtk."
#     echo "${path}"
#     echo "${data_type}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_1" -t "ttk_1" -v "${data_type}_1" -a "${anchor_pos}" -s "${figure_size}" -r "${block_1_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_2" -t "ttk_2" -v "${data_type}_2" -a "${anchor_pos}" -s "${figure_size}" -r "${block_2_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_3" -t "ttk_3" -v "${data_type}_3" -a "${anchor_pos}" -s "${figure_size}" -r "${block_3_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_1_up" -t "ttk_1_up" -v "${data_type}_1_up" -a "${anchor_pos}" -s "${figure_size}" -r "${block_1_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_2_up" -t "ttk_2_up" -v "${data_type}_2_up" -a "${anchor_pos}" -s "${figure_size}" -r "${block_2_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_3_up" -t "ttk_3_up" -v "${data_type}_3_up" -a "${anchor_pos}" -s "${figure_size}" -r "${block_3_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_4_up" -t "ttk_4_up" -v "${data_type}_4_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_5_up" -t "ttk_5_up" -v "${data_type}_5_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_6_up" -t "ttk_6_up" -v "${data_type}_6_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_7_up" -t "ttk_7_up" -v "${data_type}_7_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_8_up" -t "ttk_8_up" -v "${data_type}_8_up" -a "${anchor_pos}" -s "${figure_size}"
    
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}-1-local" -t "ttk-1-local" -a "${anchor_pos}" -s "${figure_size}" -r "${block_1_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}-2-local" -t "ttk-2-local" -a "${anchor_pos}" -s "${figure_size}" -r "${block_2_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}-3-local" -t "ttk-3-local" -a "${anchor_pos}" -s "${figure_size}" -r "${block_3_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}-1-up-local" -t "ttk-1-up-local" -a "${anchor_pos}" -s "${figure_size}" -r "${block_1_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}-2-up-local" -t "ttk-2-up-local" -a "${anchor_pos}" -s "${figure_size}" -r "${block_2_local}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}-3-up-local" -t "ttk-3-up-local" -a "${anchor_pos}" -s "${figure_size}" -r "${block_3_local}"
    


# else
#   echo "plot csv"
#   python ./python/plot_cpt.py -p "${path}" -f "${data_type}_1" -t "ttk_1" -o "${data_type}_1" -a "${anchor_pos}" -s "${figure_size}"
#   python ./python/plot_cpt.py -p "${path}" -f "${data_type}_2" -t "ttk_2" -o "${data_type}_2" -a "${anchor_pos}" -s "${figure_size}"
#   python ./python/plot_cpt.py -p "${path}" -f "${data_type}_3" -t "ttk_3" -o "${data_type}_3" -a "${anchor_pos}" -s "${figure_size}"
#   python ./python/plot_cpt.py -p "${path}" -f "${data_type}_1_up" -t "ttk_1_up" -o "${data_type}_1_up" -a "${anchor_pos}" -s "${figure_size}"
#   python ./python/plot_cpt.py -p "${path}" -f "${data_type}_2_up" -t "ttk_2_up" -o "${data_type}_2_up" -a "${anchor_pos}" -s "${figure_size}"
#   python ./python/plot_cpt.py -p "${path}" -f "${data_type}_3_up" -t "ttk_3_up" -o "${data_type}_3_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_4_up" -t "ttk_4_up" -o "${data_type}_4_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_5_up" -t "ttk_5_up" -o "${data_type}_5_up" -a "${anchor_pos}" -s "${figure_size}"
    # python ./python/plot_cpt.py -p "${path}" -f "${data_type}_9_up" -t "ttk_9_up" -o "${data_type}_9_up" -a "${anchor_pos}" -s "${figure_size}"
# fi



# python ./python/find_overlap.py -p "${path}" -f "${data_type}_1" -t "ttk_1" -o "${data_type}_1_overlap" -d "1" -g "${extension}" -o1c "${data_type}_1_with_neighbor" -o2c "ttk_1_with_neighbor" -o1f "${data_type}_1_without_neighbor" -o2f "ttk_1_without_neighbor"
# python ./python/find_overlap.py -p "${path}" -f "${data_type}_2" -t "ttk_2" -o "${data_type}_2_overlap" -d "1" -g "${extension}" -o1c "${data_type}_2_with_neighbor" -o2c "ttk_2_with_neighbor" -o1f "${data_type}_2_without_neighbor" -o2f "ttk_2_without_neighbor"
# python ./python/find_overlap.py -p "${path}" -f "${data_type}_3" -t "ttk_3" -o "${data_type}_3_overlap" -d "1" -g "${extension}" -o1c "${data_type}_3_with_neighbor" -o2c "ttk_3_with_neighbor" -o1f "${data_type}_3_without_neighbor" -o2f "ttk_3_without_neighbor"

python ./python/find_overlap.py -p "${path}" -f "${data_type}_1_up" -t "ttk_1_up" -o "${data_type}_1_overlap" -d "0.1" -g "${extension}"  -o1c "${data_type}_1_up_with_neighbor" -o2c "ttk_1_up_with_neighbor" -o1f "${data_type}_1_up_without_neighbor" -o2f "ttk_1_up_without_neighbor"
python ./python/find_overlap.py -p "${path}" -f "${data_type}_2_up" -t "ttk_2_up" -o "${data_type}_2_overlap" -d "0.1" -g "${extension}" -o1c "${data_type}_2_up_with_neighbor" -o2c "ttk_2_up_with_neighbor" -o1f "${data_type}_2_up_without_neighbor" -o2f "ttk_2_up_without_neighbor"
python ./python/find_overlap.py -p "${path}" -f "${data_type}_3_up" -t "ttk_3_up" -o "${data_type}_3_overlap" -d "0.1" -g "${extension}" -o1c "${data_type}_3_up_with_neighbor" -o2c "ttk_3_up_with_neighbor" -o1f "${data_type}_3_up_without_neighbor" -o2f "ttk_3_up_without_neighbor"


# python ./python/find_overlap.py -p "${path}" -f "${data_type}_10" -t "ttk_10" -o "${data_type}_10_overlap" -d "1" -g "${extension}" -o1c "${data_type}_10_with_neighbor" -o2c "ttk_10_with_neighbor" -o1f "${data_type}_10_without_neighbor" -o2f "ttk_10_without_neighbor"
# python ./python/find_overlap.py -p "${path}" -f "${data_type}_11" -t "ttk_11" -o "${data_type}_11_overlap" -d "1" -g "${extension}" -o1c "${data_type}_11_with_neighbor" -o2c "ttk_11_with_neighbor" -o1f "${data_type}_11_without_neighbor" -o2f "ttk_11_without_neighbor"

python ./python/find_overlap.py -p "${path}" -f "${data_type}_10_up" -t "ttk_10_up" -o "${data_type}_10_overlap" -d "0.1" -g "${extension}" -o1c "${data_type}_10_up_with_neighbor" -o2c "ttk_10_up_with_neighbor" -o1f "${data_type}_10_up_without_neighbor" -o2f "ttk_10_up_without_neighbor"
python ./python/find_overlap.py -p "${path}" -f "${data_type}_11_up" -t "ttk_11_up" -o "${data_type}_11_overlap" -d "0.1" -g "${extension}" -o1c "${data_type}_11_up_with_neighbor" -o2c "ttk_11_up_with_neighbor" -o1f "${data_type}_11_up_without_neighbor" -o2f "ttk_11_up_without_neighbor"


# python ./python/find_overlap.py -p "${path}" -f "${data_type}_4_up" -t "ttk_4_up" -o "${data_type}_4_overlap" -d "0.01" -g "${extension}"  -o1c "${data_type}_4_up_with_neighbor" -o2c "ttk_4_up_with_neighbor" -o1f "${data_type}_4_up_without_neighbor" -o2f "ttk_4_up_without_neighbor"

# python ./python/find_overlap.py -p "${path}" -f "${data_type}_5_up" -t "ttk_5_up" -o "${data_type}_5_overlap" -d "0.01" -g "${extension}"  -o1c "${data_type}_5_up_with_neighbor" -o2c "ttk_5_up_with_neighbor" -o1f "${data_type}_5_up_without_neighbor" -o2f "ttk_5_up_without_neighbor"

# python ./python/find_overlap.py -p "${path}" -f "${data_type}_6_up" -t "ttk_6_up" -o "${data_type}_6_overlap" -d "0.01" -g "${extension}"  -o1c "${data_type}_6_up_with_neighbor" -o2c "ttk_6_up_with_neighbor" -o1f "${data_type}_6_up_without_neighbor" -o2f "ttk_6_up_without_neighbor"

# python ./python/find_overlap.py -p "${path}" -f "${data_type}_12" -t "ttk_12" -o "${data_type}_12_overlap" -d "1"
# python ./python/find_overlap.py -p "${path}" -f "${data_type}_12_up" -t "ttk_12_up" -o "${data_type}_12_up_overlap" -d "0.1"



# python ./python/find_overlap.py -p "${path}" -f "${data_type}" -t "ttk" -o "${data_type}_overlap" -d "1" -g "${extension}"

# python ./python/find_overlap.py -p "${path}" -f "${data_type}" -t "ttk" -o "${data_type}_overlap" -d "1" -g "${extension}"

# python ./python/extract_mesh.py -n "${data_type}" -p "${path}"
# python ./python/create_cube.py
# echo "End Runing Script"
# # ```


# # python to create critical point
# source ~/enter/etc/profile.d/conda.sh
# conda activate mfa
# python ./python/extract_coast_line.py -p "${path}" -i "coastLine.vtp" -o "coastLine2.vtk" -r "${block_2}"
# python ./python/extract_coast_line.py -p "${path}" -i "coastLine.vtp" -o "coastLine3.vtk" -r "${block_3}"
# python ./python/extract_coast_line.py -p "${path}" -i "coastLine.vtp" -o "coastLine1.vtk" -r "${block_1}"