import matplotlib.pyplot as plt

# Data extracted from the table
folder = "../../build/examples/vis-curves/"

fig_size_s = (7, 5)
fig_size_gamma = (5, 5)
fig_size_epsilon = (5, 5)
plt.rcParams.update({'font.size': 31})
plt.rcParams["text.usetex"] = True


def draw_function(ax,xlabel, xlabel_title, num_cc, num_loops, # fig_name, pdf_name, 
                  legend_loc, bbox_to_anchor, fig_size,x_label_rotation=0): 
    x_axes = range(len(xlabel))

    # fig, ax = plt.subplots(figsize=fig_size)


    
    ax.plot(num_loops, marker='o', markersize=15, label=r"\#Loop")
    ax.plot(num_cc, marker='s', markersize=15, label=r"\#CC")

    ax.set_ylim(-0.1*max(max(num_loops), max(num_cc)), 1.1*max(max(num_loops), max(num_cc)))
    ax.set_xticks(x_axes)
    # ax.set_xticklabels(xlabel)
    
    ax.set_xticklabels(xlabel, rotation=x_label_rotation)
    
    ax.set_xlabel(xlabel_title)

    ax.legend(loc=legend_loc, bbox_to_anchor=bbox_to_anchor)
    ax.grid(False)

    plt.setp(ax.get_lines(), linewidth=6)
    plt.setp(ax.get_legend().get_lines(), linewidth=6)


    # Use tight_layout to preserve labels and axes
    # fig.tight_layout(pad=0.1)


    # # Save without white border
    # fig.savefig(pdf_name, bbox_inches='tight', pad_inches=0)
    # fig.savefig(fig_name, bbox_inches='tight', pad_inches=0)


    # plt.close(fig)

# def draw_function(xlabel, xlabel_title, num_cc, num_loops, fig_name, pdf_name, legned_loc,bbox_to_anchor, fig_size):
#     x_axes = range(len(xlabel))
#     # Plotting
#     plt.figure(figsize=fig_size)

#     # plt.rcParams.update({'font.size': 25})

#     plt.plot(num_loops, marker='o', markersize=15, label='#Loop')
#     plt.plot(num_cc, marker='s',markersize=15, label='#CC')

#     plt.ylim(0, max(max(num_loops), max(num_cc))+5)

#     plt.xticks(x_axes, xlabel)  # Update the x-axis labels

#     # plt.xlabel(r"step size $s$")
#     plt.xlabel(xlabel_title)

#     # plt.ylabel('Count')
#     # plt.title('Number of Loops and Connected Components vs. $\gamma$')
#     plt.legend(loc=legned_loc, bbox_to_anchor=bbox_to_anchor)
#     # plt.legend(loc='lower right', bbox_to_anchor=(1, 0.1))

#     # plt.axvline(x=chosen_pos, linestyle='--', color='red', label=r'Selected $s$')
#     # plt.text(dash_line_pos, 0, r"Selected $s$", ha=dash_line_loc, va='bottom', color='red')
#     plt.grid(False)
#     plt.tight_layout()
#     plt.setp(plt.gca().get_lines(), linewidth=6)
#     plt.setp(plt.gca().get_legend().get_lines(), linewidth=6)
#     plt.savefig(pdf_name)
#     plt.savefig(fig_name)


# Data for the plots
folder = "../../build/examples/vis-curves/"


def sinc_fig():

    # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]
    num_cc_contour = [32, 32, 32, 32]
    num_loops_contour = [20, 28, 28, 28]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [32, 32, 32, 32]
    num_loops_epsilon = [28, 28, 28, 28]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [32, 32, 32, 32]
    num_loops_gamma = [24, 28, 28, 28]
    xlabel_title_gamma = r"Connection threshold $\gamma$"


    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 subplots in a row

    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.1), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'sinc_contour_para.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'sinc_contour_para.pdf', bbox_inches='tight', pad_inches=0)



def sinc_1_fig():

    # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]
    num_cc_contour = [67,60,60,60]
    num_loops_contour = [37,40,40,40]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [60, 60, 60, 60]
    num_loops_epsilon = [40, 40, 40, 40]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [65, 60, 60, 60]
    num_loops_gamma = [33,39,40,40]
    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 subplots in a row

    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.1), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_sinc_contour_0.33.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_sinc_contour_0.33.pdf', bbox_inches='tight', pad_inches=0)

def s3d_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour = [21, 24, 23, 22, 22, 22, 22,22]
    num_loops_contour = [18,19,17,18,18,18, 18,18]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [18, 18, 18, 18]
    num_loops_epsilon = [22,22,22,22]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [17, 18, 18, 18]
    num_loops_gamma = [25,22,22,22]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.1), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 's3d_contour_para.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 's3d_contour_para.pdf', bbox_inches='tight', pad_inches=0)
    

def s3d_1_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour = [23,24,23,23,23,23,23,23]
    num_loops_contour = [22,22,21,21,21,21,21,21]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [23, 23, 23, 23]
    num_loops_epsilon = [21,21,21,21]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [24,23,23,23]
    num_loops_gamma = [19,21,21,21]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.1), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_s3d_contour_30.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_s3d_contour_30.pdf', bbox_inches='tight', pad_inches=0)
       


def s3d_3_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour = [21,21,21,21,21,21,21,21]
    num_loops_contour = [13,13,13,13,13,13,13,13]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [21,21,21,21]
    num_loops_epsilon = [13,13,13,13]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [25,21,21,21]
    num_loops_gamma = [12,13,13,13]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.025), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.025), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.025), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_s3d_contour_60.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_s3d_contour_60.pdf', bbox_inches='tight', pad_inches=0)
       


def schwefel_1_fig():
    # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]
    num_cc_contour =  [39, 39, 39, 39]
    num_loops_contour = [32, 32, 32, 32]
    xlabel_title_contour = r"Step size $s$"


    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [39, 39, 39, 67]
    num_loops_epsilon = [32, 32, 32, 47]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [54,39, 39, 39]
    num_loops_gamma = [20,31,32,32]
    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 subplots in a row

    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (0.8, 0.025), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.025), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_schwefel_contour_100.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_schwefel_contour_100.pdf', bbox_inches='tight', pad_inches=0)    

def schwefel_2_fig():
    # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]
    num_cc_contour =  [21, 22, 22, 22]
    num_loops_contour = [19, 21, 21, 21]
    xlabel_title_contour = r"Step size $s$"


    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [22, 22, 22, 75]
    num_loops_epsilon = [21, 21, 21, 54]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [52, 24, 22, 22]
    num_loops_gamma = [4, 17, 21, 21]
    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 subplots in a row

    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (0.8, 0.6), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1.05, 0.6), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_schwefel_contour_500.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_schwefel_contour_500.pdf', bbox_inches='tight', pad_inches=0)

def gaussian_pair_fig():
    # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]
    num_loops_contour =  [1, 0, 0, 0]
    num_cc_contour = [3,2,2,2]
    xlabel_title_contour = r"Step size $s$"


    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [7,4,2,2]
    num_loops_epsilon = [1,0,0,0]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [2,2,2,2]
    num_loops_gamma = [0,0,0,0]
    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 subplots in a row

    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.65), fig_size_epsilon)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.65), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.3), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_gaussian_pair_js.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_gaussian_pair_js.pdf', bbox_inches='tight', pad_inches=0)


def karman_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour = [106, 133, 129, 127, 124, 124, 124, 124]
    num_loops_contour = [81, 89, 85, 82, 82, 82, 82, 82]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [138, 124, 124, 124]
    num_loops_epsilon =  [97, 82, 82, 82]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [124, 124, 124, 124]
    num_loops_gamma = [82, 82, 82, 82]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.05), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.05), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.05), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_karman.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_karman.pdf', bbox_inches='tight', pad_inches=0)

def hurricane_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour =[1920, 2067, 2067, 2038, 2012, 2015, 2015, 2015]
    num_loops_contour =  [1743, 1850, 1926, 1914, 1905, 1901, 1901, 1901]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon =[2014, 2015, 2015, 2015]
    num_loops_epsilon =  [1903, 1901, 1901, 1901]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [2022,2015, 2015, 2014]
    num_loops_gamma = [1900, 1900, 1901, 1902]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.1), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.1), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_hurricane.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_hurricane.pdf', bbox_inches='tight', pad_inches=0)


def gaussian_mix_fig():
    # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]
    num_loops_contour =  [12, 2,2, 2]
    num_cc_contour = [2,1,1,1]
    xlabel_title_contour = r"Step size $s$"


    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [3,1,1,1]
    num_loops_epsilon = [4,3,2,2]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [7,4,1,1]
    num_loops_gamma = [1,2,2,7]
    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 subplots in a row

    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.65), fig_size_epsilon)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1.08, 0.68), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (0.9, 0.68), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_gaussian_mixture.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_gaussian_mixture.pdf', bbox_inches='tight', pad_inches=0)


def cesm_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour =[788, 910, 948, 912, 860, 841, 851, 869]
    num_loops_contour = [4713, 4245, 3679, 3477, 3402, 3334, 3250, 3283]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon =[58, 57, 57, 57]
    num_loops_epsilon =  [31, 31, 31, 31]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma = [60, 57, 57, 57]
    num_loops_gamma = [31, 31, 31, 31]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.195), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.02), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.02), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_cesm.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_cesm.pdf', bbox_inches='tight', pad_inches=0)


def boussinesq_fig():
     # Contour data
    xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    num_cc_contour =[493, 483, 493, 485, 479, 479, 479, 479]
    num_loops_contour = [548, 452, 401, 405, 398, 398, 398, 398]
    xlabel_title_contour = r"Step size $s$"

    # Epsilon data
    xlabel_epsilon = [r'1$e^{-6}$', r'1$e^{-8}$', r'1$e^{-10}$', r'1$e^{-12}$']
    num_cc_epsilon = [1545, 479, 479, 498]
    num_loops_epsilon =  [959, 398, 398, 417]
    xlabel_title_epsilon = r"Accuracy threshold $\epsilon$"

    # Gamma data
    xlabel_gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']
    num_cc_gamma =  [483, 480, 479, 479]
    num_loops_gamma =  [392, 395, 398, 398]

    xlabel_title_gamma = r"Connection threshold $\gamma$"

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(17, 5),gridspec_kw={'width_ratios': [7, 5, 5]})  # 3 subplots in a row
 # Relative widths of the subplots


    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, num_cc_contour, num_loops_contour, 'lower right', (1, 0.1), fig_size_s,30)
    draw_function(axes[1], xlabel_epsilon, xlabel_title_epsilon, num_cc_epsilon, num_loops_epsilon, 'lower right', (1, 0.5), fig_size_epsilon)
    draw_function(axes[2], xlabel_gamma, xlabel_title_gamma, num_cc_gamma, num_loops_gamma, 'lower right', (1, 0.1), fig_size_gamma)

    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'para_boussinesq.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'para_boussinesq.pdf', bbox_inches='tight', pad_inches=0)




s3d_fig()
s3d_fig()
s3d_3_fig()
s3d_1_fig()
sinc_fig()
sinc_1_fig()
schwefel_1_fig()
schwefel_2_fig()
gaussian_pair_fig()
karman_fig()
hurricane_fig()
gaussian_mix_fig()
cesm_fig()
boussinesq_fig()



# sinc_contour_epsilon()
# # sinc_contour_epsilon()
# sinc_contour_s()

# # sinc_contour_gamma()

# s3d_contour_s()
# s3d_contour_epsilon()
# s3d_contour_gamma()



# #sinc - contour
#step size

# #step size
# step_size = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16"]



# ##epsilon
# epsilon = 
# gamma = [r'1.0$s$', r'1.5$s$', r'2.0$s$', r'2.5$s$']

#     draw_function(xlabel, r"step size $s$", num_cc, num_loops, fig_name, pdf_name, 'upper right', (1, 0.9))




#epsilon





#S3D - contour
# num_loops = [18, 19, 17, 18, 18, 18]
# gamma = [r'1.0$\gamma$', r'1.5$\gamma$', r'2.0$\gamma$', r'2.5$\gamma$']
# num_cc = [21,24,23,22,22,22]
# x_axes = [0, 1, 2, 3, 4, 5]
# step_size = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64"]
# fig_name=folder + 'S3D-contour-step-size.png'
# pdf_name=folder + 'S3D-contour-step-size.pdf'
# chosen_pos = 3
# dash_line_loc = 'right'
# dash_line_pos=chosen_pos-0.05


