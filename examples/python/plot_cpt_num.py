import matplotlib.pyplot as plt

x = [1, 2, 3, 4, 5, 6,7,8,9]
x_label = ['1x', '$10^2$x', '$20^2$x', '$30^2$x', '$40^2$x', '$50^2$x', '$60^2$x', '$70^2$x', '$80^2$x']

MFA_1=[114,123,124,124,124,124,124,124,124]
ttk_1=[107,134,136,150,166,188,226,283,381]
# ttk_3=[74, 79, 79, 79, 79, 80]
original_1=[53,53,53,53,53,53,53,53,53]
# original_2=[53,53,53,53,53,53]
# original_3=[43,43,43,43,43,43]

# MFA_1=[124,124,124,124,124,124]
# MFA_2=[116,116,116,116,116,116]
# MFA_3=[78,78,78,78,78,78]



plt.figure(figsize=(9, 5))

# plt.plot(x,ttk_3, color = (1.0, 0.647, 0.0, 1.0), marker='x', label = 'TTK Threshold 1e-3')
# # plt.plot(x,ttk_2, color = 'red', marker='o', label = 'Threshold 2e-4')
# plt.plot(x,original_3, color = (0.0, 0.5, 0.0, 1.0), marker='o', label = 'Original Threshold 1e-3')
# plt.plot(x,MFA_3, color = (0.0, 0.0, 1.0, 1.0), marker='^', label = 'MFA Threshold 1e-3')

# plt.plot(x,ttk_2, color = (1.0, 0.647, 0.0, 0.5), marker='x', label = 'Threshold 2e-4')
# plt.plot(x,original_2, color = (0.0, 0.5, 0.0, 0.5), marker='o', label = 'Original Threshold 2e-4')
# plt.plot(x,MFA_2, color =(0.0, 0.0, 1.0, 0.5), marker='^', label = 'MFA Threshold 2e-4')

plt.plot(x,ttk_1, color = (1.0, 0.647, 0.0, 1.0), marker='x', label = 'MFA-TTK')
plt.plot(x,original_1, color = (0.0, 0.5, 0.0, 1.0), marker='o', label = 'raw data-TTK')
plt.plot(x,MFA_1, color = (0.0, 0.0, 1.0, 1.0), marker='^', label = 'MFA-our method')




# Adjust layout to make room for the legend


plt.xticks(x, x_label)
plt.xlabel('Upsampling Ratio')
plt.ylabel('Critical Point Number')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.ylim(bottom=40, top = 400)

plt.savefig('cpt_comp.png', dpi=600) # Saves the figure with high resolution
# plt.show()