import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
import pyvista as pv
import argparse
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from itertools import product
import matplotlib.lines as mlines
from PIL import Image
import os
from matplotlib.patches import Polygon



def draw_quad(ax,center,length):
    p1=[center[0]-length/2,center[1]-length/2,center[2]]
    p2=[center[0]+length/2,center[1]-length/2,center[2]]
    p3=[center[0]+length/2,center[1]+length/2,center[2]]
    p4=[center[0]-length/2,center[1]+length/2,center[2]]
    
    edges = [[p1, p2, p3, p4, p1]]
    ax.add_collection3d(Line3DCollection(edges, colors=(0.5,0.0,0.8), linewidths=1, linestyles='-'))

def draw_line(ax,center,length):
    p1=[center[0]-length/2,center[1]-length/2,center[2]]
    p2=[center[0]+length/2,center[1]+length/2,center[2]]
    
    p3=[center[0]-length/2,center[1]+length/2,center[2]]
    p4=[center[0]+length/2,center[1]-length/2,center[2]]
    
    edges = [[p1, p2],[p3, p4]]
    ax.add_collection3d(Line3DCollection(edges, colors='black', linewidths=1, linestyles='-'))
    
    

def draw_cube(ax,center, size):
    # Generate the vertices of the cube
    r = [-size/2, size/2]
    vertices = np.array(list(product(r, r, r))) + center
    # Generate the sides of the cube
    edges = [[vertices[i] for i in [0,1,3,2,0,4,5,7,6,4,5,1,3,7,6,2]]]
    ax.add_collection3d(Line3DCollection(edges, colors='blue', linewidths=0.5, linestyles='-'))


def draw_a_square(ax,center, length_x, length_y, length_z,color):
    vertices = np.array([
    [center[0] + length_x / 2, center[1] - length_y / 2, center[2] - length_z / 2],
    [center[0] + length_x / 2, center[1] + length_y / 2, center[2] - length_z / 2],
    [center[0] + length_x / 2, center[1] + length_y / 2, center[2] + length_z / 2],
    [center[0] + length_x / 2, center[1] - length_y / 2, center[2] + length_z / 2]
    ])
    edges=[[vertices[0],vertices[1],vertices[2],vertices[3],vertices[0]]]
    
    ax.add_collection3d(Line3DCollection(edges, colors=color, linewidths=1.5,linestyles='-'))
    


def draw_front_cube_different_length(ax,center,length_x, length_y, length_z,color):
    vertices = np.array([
    [center[0] - length_x / 2, center[1] - length_y / 2, center[2] - length_z / 2],
    [center[0] + length_x / 2, center[1] - length_y / 2, center[2] - length_z / 2],
    [center[0] + length_x / 2, center[1] + length_y / 2, center[2] - length_z / 2],
    [center[0] - length_x / 2, center[1] + length_y / 2, center[2] - length_z / 2],
    [center[0] - length_x / 2, center[1] - length_y / 2, center[2] + length_z / 2],
    [center[0] + length_x / 2, center[1] - length_y / 2, center[2] + length_z / 2],
    [center[0] + length_x / 2, center[1] + length_y / 2, center[2] + length_z / 2],
    [center[0] - length_x / 2, center[1] + length_y / 2, center[2] + length_z / 2]
    ])
    
    # edges = [
    # [vertices[0], vertices[1], vertices[2], vertices[3], vertices[0]],  # Bottom
    # [vertices[4], vertices[5], vertices[6], vertices[7], vertices[4]],  # Top
    # [vertices[0], vertices[4]], [vertices[1], vertices[5]],  # Sides
    # [vertices[2], vertices[6]], [vertices[3], vertices[7]]
    # ]
    # edges=[[vertices[1],vertices[5],vertices[6],vertices[2],vertices[1]],
    #        [vertices[6],vertices[7],vertices[3],vertices[2]],
    #        [vertices[5],vertices[4],vertices[7]]]
    edges=[[vertices[1],vertices[5],vertices[4],vertices[0],vertices[1]],
           [vertices[1],vertices[2],vertices[6],vertices[5]],
           [vertices[4],vertices[7],vertices[6]]]
    
    ax.add_collection3d(Line3DCollection(edges, colors=color, linewidths=1.5,linestyles='-'))
    
def draw_back_cube_different_length(ax,center,length_x, length_y, length_z, color):
    vertices = np.array([
    [center[0] - length_x / 2, center[1] - length_y / 2, center[2] - length_z / 2],
    [center[0] + length_x / 2, center[1] - length_y / 2, center[2] - length_z / 2],
    [center[0] + length_x / 2, center[1] + length_y / 2, center[2] - length_z / 2],
    [center[0] - length_x / 2, center[1] + length_y / 2, center[2] - length_z / 2],
    [center[0] - length_x / 2, center[1] - length_y / 2, center[2] + length_z / 2],
    [center[0] + length_x / 2, center[1] - length_y / 2, center[2] + length_z / 2],
    [center[0] + length_x / 2, center[1] + length_y / 2, center[2] + length_z / 2],
    [center[0] - length_x / 2, center[1] + length_y / 2, center[2] + length_z / 2]
    ])
    
    # edges = [
    # [vertices[0], vertices[1], vertices[2], vertices[3], vertices[0]],  # Bottom
    # [vertices[4], vertices[5], vertices[6], vertices[7], vertices[4]],  # Top
    # [vertices[0], vertices[4]], [vertices[1], vertices[5]],  # Sides
    # [vertices[2], vertices[6]], [vertices[3], vertices[7]]
    # ]
    edges=[[vertices[7],vertices[3]],
           [vertices[0],vertices[3],vertices[2]]]
    
    
    ax.add_collection3d(Line3DCollection(edges, colors=color, linewidths=1.5,linestyles='-'))

def plot_3d(data_name, ttk_name, path,vtk_name,fig_size,bbox_to_anchor, local_range):
    figure_name = path+data_name+'.png'
    MFA_file= path+data_name+'.csv'
    ttk_file= path+ttk_name+'.csv'       
    
    range =  tuple(map(float, local_range.split('-')))
    if data_name[:3]=="qmc":
        data_size=[68,68,114]
    elif data_name[:3]=="rti":
        data_size=[143,255,255]
    min_x=data_size[0]*range[0]
    min_y=data_size[1]*range[2]
    min_z=data_size[2]*range[4]
    max_x=data_size[0]*range[1]
    max_y=data_size[1]*range[3]
    max_z=data_size[2]*range[5]
    
    bounds=[min_x,max_x,min_y,max_y,min_z,max_z]
    center_x= (bounds[0]+bounds[1])/2
    center_y= (bounds[2]+bounds[3])/2
    center_z= (bounds[4]+bounds[5])/2    
    center_local=[center_x,center_y,center_z]
    length_local=[bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]]
    
    if vtk_name=='default':

        bounds=[min_x,max_x,min_y,max_y,min_z,max_z]
    else:
        vtk_name=path+vtk_name+'.vtk'
        volume = pv.read(vtk_name)
        bounds = volume.bounds   
        
    length = 1
    
    if data_name[-5:]== "local":
        if data_name[:3]=="qmc":
            length = 0.4
        else:
            length = 0.2
    
    if data_name=="qmcpack_4_up" or data_name=="qmcpack_5_up":
        length = 0.5
    if data_name=="rti_4_up" or data_name=="rti_5_up" or data_name=="rti_6_up" or data_name=="rti_7_up" or data_name=="rti_8_up":
        length = 0.5
    
    ttk_use=['PositionX','PositionY','PositionZ']
    MFA_use=['x0','x1','x2']
    
    critical_points_ttk = pd.read_csv(ttk_file,usecols=ttk_use).to_numpy()
    critical_points_MFA = pd.read_csv(MFA_file,usecols=MFA_use).to_numpy()
    
    center_x= (bounds[0]+bounds[1])/2
    center_y= (bounds[2]+bounds[3])/2
    center_z= (bounds[4]+bounds[5])/2
    
    
    center=[center_x,center_y,center_z]
    
    

    fig = plt.figure(figsize=fig_size)
    
    print(fig_size)
    
    ax = fig.add_subplot(111, projection='3d')
    
    draw_block=False
    
    if draw_block:
        if data_name=="rti_2_up":
            # center_1=[58,24,161]
            center_2=[73,27,157]    
            draw_front_cube_different_length(ax, center_2, 2,2,2,'blue')
            # draw_a_square(ax,center_2, 2, 2, 2,'blue')
            draw_back_cube_different_length(ax, center_2, 2,2,2,'blue')
            
        if data_name=="rti_3_up":
            center_1=[60,62,205]
            draw_front_cube_different_length(ax, center_1, 2,2,2,'blue')  
              
    
     #draw boundary
    draw_front_cube_different_length(ax, center, bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4],'black')
   
    if draw_block:
        if data_name[:3]=="qmc":
            if int(data_name[8])<4:
                draw_front_cube_different_length(ax, center_local,length_local[0],length_local[1],length_local[2] ,'red')  
        if data_name[:3]=="rti":
            if int(data_name[4])<4:
                draw_front_cube_different_length(ax, center_local,length_local[0],length_local[1],length_local[2] ,'red')  
  
        if data_name=="qmcpack_2_up":
            center_c=[37,33,76]
            draw_front_cube_different_length(ax, center_c, 2, 2, 2,'blue')    
            
            
        if data_name=="qmcpack_3_up":
            center_c=[61,20,58]
            draw_front_cube_different_length(ax, center_c, 2, 2, 2,'blue')    
        if data_name=="rti_1_up":
            # center_1=[85,66,89]
            center_2=[80,66,89.3]
            # 79, 65, 88.3
            # center_3=[72.8,68,92.8]
            draw_front_cube_different_length(ax, center_2, 2,2,2,'blue')
            # draw_front_cube_different_length(ax, center_2, 2,2,2,'blue')



    # if data_name=="rti_3_up":
    #     center_1=[51.4,51.1,204.1]
    #     draw_front_cube_different_length(ax, center_1, 0.5,0.5,0.5,'green')
        
    #
    for point in critical_points_ttk:
        # draw_cube(ax, point, length)   
        draw_quad(ax,point,length)

    for point in critical_points_MFA:
        draw_line(ax,point,length)
        
        
        
    draw_back_cube_different_length(ax, center, bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4],'black')
    
    if draw_block:
        if data_name[:3]=="qmc":
            if int(data_name[8])<4:
                draw_back_cube_different_length(ax, center_local,length_local[0],length_local[1],length_local[2] ,'red')  
        if data_name[:3]=="rti":
            if int(data_name[4])<4:
                draw_back_cube_different_length(ax, center_local,length_local[0],length_local[1],length_local[2] ,'red')  
        if data_name=="qmcpack_2_up" or data_name=="qmcpack_3_up":
            draw_back_cube_different_length(ax, center_c, 2, 2, 2,'blue')
        if data_name=="rti_1_up":
            draw_back_cube_different_length(ax, center_2, 2,2,2,'blue')
            # draw_back_cube_different_length(ax, center_2, 2,2,2,'blue')

        # if data_name=="rti_2_up":
        #     draw_back_cube_different_length(ax, center_1, 2,2,2,'blue')
    #draw boundary
    
    
    

    max_range = np.array([bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]]).max() / 2.0
    

    ax.set_xlim([center_x -max_range, center_x+max_range])
    ax.set_ylim([center_y -max_range, center_y+max_range])
    ax.set_zlim([center_z -max_range, center_z+max_range])
    
    ax.set_box_aspect([1,1,1])
    # ax.view_init(elev=90)#, azim=azimuth_angle
    
    plt.axis('off')
    cube_legend = mlines.Line2D([], [], color=(0.5,0.0,0.8), marker='s', linestyle='None',markerfacecolor='none',
                            markersize=10, label='TTK-MFA')
    line_legend = mlines.Line2D([], [], color='black', marker='x', linestyle='None',
                            markersize=10, label='CPE-MFA')

    
    # plt.legend(loc='upper right',bbox_to_anchor=bbox_to_anchor, fontsize=14)
    # Adding a color bar to indicate the values of the z dimension
    # plt.gca().set_facecolor('none')
    # plt.axis('equal')
    # plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    # Show the plot
    # plt.show()
    # if data_name[-5:]!="local":
    #     ax.legend(handles=[cube_legend, line_legend],loc=bbox_to_anchor, fontsize=14)#,
    
    plt.savefig(figure_name, dpi=300,transparent=True) # Saves the figure with high resolution
    
    

def plot_3d_point(data_name, ttk_name, path,vtk_name):
    # path = '../../build/examples/qmcpack/'    
    # folder='threshold1e_3/'
    
    figure_name = path+data_name+'.png'
    MFA_file= path+data_name+'.csv'
    ttk_file= path+ttk_name+'.csv'       
    vtk_name=path+vtk_name+'.vtk'
    
    # vtk_file_path=path+'qmcpack_'+str(index)+'.vtk'   
    volume = pv.read(vtk_name)
    
    print("run here")
    
    ttk_use=['PositionX','PositionY','PositionZ']
    MFA_use=['x0','x1','x2']
    
    critical_points_ttk = pd.read_csv(ttk_file,usecols=ttk_use).to_numpy()
    critical_points_MFA = pd.read_csv(MFA_file,usecols=MFA_use).to_numpy()
    
    
    bounds = volume.bounds

    cube = pv.Box(bounds=bounds)


    plotter = pv.Plotter(off_screen=True,window_size=(1920, 1080))
    cube_size = 0.5  # Adjust based on your preference
    tetra_size = 0.5  # Adjust based on your preference
    line_length = 2.0

    opacity_transfer_function = [(0, 0.0),  # At the lowest value, opacity is 0 (fully transparent)
                             (0.25, 0.1),  # Increase opacity slightly
                             (0.5, 0.3),  # Mid values have higher opacity
                             (0.75, 0.1),  # Decrease opacity again
                             (1, 0.0)]  # At the highest value, opacity is back to 0

    # plotter.add_volume(volume, opacity='foreground', cmap="coolwarm")
    plotter.add_mesh(cube, color='lightblue', style='wireframe')
    # Draw points from file A as cubes
    
    
    line_points = []
    lines = []

    idx=0
    for point in critical_points_MFA:
        start_x = [point[0] - line_length / 2, point[1], point[2]]
        end_x = [point[0] + line_length / 2, point[1], point[2]]
        start_y = [point[0], point[1] - line_length / 2, point[2]]
        end_y = [point[0], point[1] + line_length / 2, point[2]]
        line_points.append(start_x)
        line_points.append(end_x)
        line_points.append(start_y)
        line_points.append(end_y)
        lines.append([2, 4 * idx, 4 * idx + 1])  # Line 1
        lines.append([2, 4 * idx + 2, 4 * idx + 3])  # Line 2
        idx+=1

    line_points = np.array(line_points)
    lines = np.hstack(lines)
    # lines = np.hstack([[2, 0, 1]])
    line_poly = pv.PolyData()
    line_poly.points = line_points
    line_poly.lines = lines

    
    # print(lines)
    # print(line_points)
    
    # lines = np.hstack(lines)
    
    print("run here 1")
    
    # line_poly= pv.PolyData(line_points,lines=lines)
    
    
    print("run here 2")    #run line 2
    
    plotter.add_mesh(line_poly, color="black", line_width=2)
    
    print("run here 3")    
    # for point in critical_points_MFA:
    #     cube = pv.Cube(center=point, x_length=cube_size, y_length=cube_size, z_length=cube_size)
    #     plotter.add_mesh(cube, color='black',style='wireframe')  # You can change the color
    
    for point in critical_points_ttk:
        tetra = pv.Tetrahedron(center=point, radius=tetra_size)
        tetra.opacity = 0.0
        plotter.add_mesh(tetra, color='red',style='wireframe',line_width=2)  

        
    print("run here 5")   
    # pv.global_theme.axes.show = False


    plotter.show()
    print(figure_name)
    plotter.screenshot(figure_name, transparent_background=True)#, return_img=False #window_size=[600, 600], 
    plotter.close()
    print("run here 4")    
    

def plot_fuc(data_name, ttk_name, path,ply_name,bbox_to_anchor,fig_size):     
    
    is_block=False
    
    # ttk_name='ttk_1'
    # data_name='cesm_1'
    
    # path = '../../build/examples/cesm/'    
    # folder='threshold1e_3/'
    
    # draw_block=1
    
    # if draw_block==1:
    #     block='_block_1'
    # elif draw_block==2:
    #     block='_block_2_'
    
    
    
    # if is_block:
    #     file_name = path + data_name+block+'.ply'
    #     figure_name = path+data_name+block+'.png'
    #     MFA_file= path+data_name+block+'.csv'
    #     ttk_file= path+ttk_name + block +'.csv'
    # else:
    file_name = path + ply_name+'.ply' #up+
    figure_name = path+data_name+'.png'
    MFA_file= path+data_name+'.csv'
    ttk_file= path+ttk_name+'.csv'        

    # if index == 50:
    #     file_name = path + 'up_'+str(index-10)+'.ply'
    # else:
    # file_name = path + 'patch_'+str(index)+'.ply'


    # path_2 = path #+folder

    # orignal_file= path_2+'original.csv'
    # simplification_file = path+'simplification_30.csv'

    plt.clf()
    plt.figure(figsize=fig_size)#
    
    
    
    
    mesh = pv.read(file_name)
    vertices = mesh.points
    faces = mesh.faces.reshape((-1, 4))[:, 1:4]

    # Extract x, y coordinates and z values
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]

    z_faces = np.mean(z[faces], axis=1)
    
    # max=vertices.max(axis=0)
    # min=vertices.min(axis=0)    
    # bbox = Rectangle((min[0], min[1]), max[0]-min[0], max[1]-min[1], fill=False, edgecolor='white', linewidth=10)
    
    # plt.gca().add_patch(bbox)
    

    if data_name[:4]=="cesm":
        vmin=-1
        vmax=430
    if data_name[:4]=="s3d_":
        vmin=-1
        vmax=360


    # Create a Triangulation object from the vertices and faces
    triangulation = tri.Triangulation(x, y, triangles=faces)
    # Draw triangles with color based on the z value
    tripcolor_plot= plt.tripcolor(triangulation, facecolors=z_faces, cmap='coolwarm', edgecolors='none',vmin=vmin, vmax=vmax,  rasterized=True) #summer # bwr




    critical_points_ttk = pd.read_csv(ttk_file)
    critical_points_MFA = pd.read_csv(MFA_file)
    
    # bbox_to_anchor=(1.0, 1.07)
    # if is_block==False:
    # bbox_to_anchor=(1.0, 1.0)
    scale = 25
    line_width=1
    
    
    is_block = False
    
    # if is_block == False:
    #     rectangle = Rectangle(p1, p2[0] - p1[0], p2[1] - p1[1], edgecolor='orange', facecolor='none')
    #     plt.gca().add_patch(rectangle)
        # if (data_name =='cesm' and index !=2) or (data_name =='s3d' and index !=1):
        #     rectangle2 = Rectangle(p3, p4[0] - p3[0], p4[1] - p3[1], edgecolor='pink', facecolor='none')
        #     plt.gca().add_patch(rectangle2)

    if is_block:  
        if data_name =="cesm_1_up":
            block_vertices = [(180,1266),(180,1284),(216,1284),(216,1266)]
            block = Polygon(block_vertices, closed=True,linewidth=3,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block)
            block_vertices2=[(277,1301),(277,1319),(292,1319),(292,1301)]
            block2 = Polygon(block_vertices2, closed=True,linewidth=3,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block2)
        
    edge_color=(0.5,0.0,0.8)
    # edge_color=(75/255,0,146/255)

    plt.scatter(critical_points_ttk['PositionX'], critical_points_ttk['PositionY'], s=scale, c='none', marker='s', label='TTK-MFA', edgecolor=edge_color,linewidths=line_width) #green (0.004,0.81,0.15) #purple (0.8,0.0,0.89)  (0.8,0.0,0.89)
    plt.scatter(critical_points_MFA['x0'], critical_points_MFA['x1'], s=scale, c='black', marker='x', label='CPE-MFA',linewidths=line_width)
    plt.axis('off')
    
    
    plt.legend(loc='upper right',bbox_to_anchor=bbox_to_anchor, fontsize=14)
    # Adding a color bar to indicate the values of the z dimension
    # plt.colorbar(label='Function Value (z)')
    if is_block:    
        line_width=1.5
        if data_name =="s3d_3_up":
            block_vertices = [(160,101),(160,107),(164,107),(164,101)]
            # block_vertices = [(158,99),(158,109),(166,109),(166,99)]
            block = Polygon(block_vertices, closed=True,linewidth=line_width,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block)
            block_vertices2=[(195,108),(195,115),(201,115),(201,108)]
            # block_vertices2=[(193,106),(193,117),(203,117),(203,106)]
            block2 = Polygon(block_vertices2, closed=True,linewidth=line_width,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block2)
            block_vertices3=[(98,57),(98,62),(102,62),(102,57)]
            # block_vertices3=[(96,55),(96,64),(104,64),(104,55)]
            block3 = Polygon(block_vertices3, closed=True,linewidth=line_width,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block3)
            block_vertices4=[(158,58),(158,64),(166,64),(166,58)]
            # block_vertices4=[(156,56),(156,66),(168,66),(168,56)]
            block4 = Polygon(block_vertices4, closed=True,linewidth=line_width,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block4)
            block_vertices5=[(216,53),(216,57),(223,57),(223,53)]
            # block_vertices5=[(214,51),(214,59),(225,59),(225,51)]
            block5 = Polygon(block_vertices5, closed=True,linewidth=line_width,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block5)
        if data_name == "s3d_2_up":
            block_vertices = [(198,203),(198,209),(202,209),(202,203)]
            # block_vertices = [(196,201),(196,211),(204,211),(204,201)]
            block = Polygon(block_vertices, closed=True,linewidth=line_width,  edgecolor='white', facecolor='none')
            plt.gca().add_patch(block)

    
    plt.gca().set_facecolor('none')
    plt.colorbar(tripcolor_plot,shrink=0.7, ax=plt.gca(), location='right')
    plt.axis('equal')
    plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    # Show the plot
    plt.savefig(figure_name, dpi=300,transparent=True) # Saves the figure with high resolution


# def plot_fuc(data_name, ttk_name, path, ply_name):     

def cropping(file_name,obj_name):
    
    if obj_name[:3].lower()=="qmc":
        if int(obj_name[8])<4:
            crop_box = (796,531, 2776,2531)  # for qmcpack
            image = Image.open(file_name)

            # Crop the image and save it
            cropped_image = image.crop(crop_box)
            cropped_image.save(file_name)
    else:
        if int(obj_name[4])<4: #rti
            crop_box = (265,104, 1745,1524)   
            image = Image.open(file_name)
            
            # Crop the image and save it
            cropped_image = image.crop(crop_box)
            cropped_image.save(file_name)


    


parser = argparse.ArgumentParser(description='plot critical points.')

parser.add_argument('-p', '--path', type=str, default='path', help='path to all files')
parser.add_argument('-f', '--data_name', type=str, default='cesm_1', help='cesm_1, cesm_2...')
parser.add_argument('-t','--ttk_name', type=str, default='ttk_1', help='ttk critical point file name')
parser.add_argument('-o','--ply_name', type=str, default='default', help='cesm_1, cesm_1_up')
parser.add_argument('-v','--vtk_name', type=str, default='default', help='vtk file name')
parser.add_argument('-a',"--bbox_to_anchor", type=str,default='1.0,1.0', help='Input legend bbox_to_anchor in tuple format')
parser.add_argument('-s',"--fig_size", type=str,default='1,1', help='Input fig_size in tuple format')
parser.add_argument('-r', '--local_range', type=str, default='0-1-0-1-0-1', help='the range of points')

args = parser.parse_args()


path=args.path
data_name=args.data_name
ttk_name=args.ttk_name
ply_name=args.ply_name
vtk_name = args.vtk_name
local_range =  args.local_range

# path = "../../build/examples/qmcpack/"
# data_name = "qmcpack_1"
# ttk_name = "ttk_1"
# ply_name = "qmcpack_1"
# vtk_name = "qmcpack_1"

bbox_to_anchor = tuple(map(float, args.bbox_to_anchor.split(',')))
fig_size = tuple(map(float, args.fig_size.split(',')))



if ply_name=='default':
    plot_3d(data_name, ttk_name, path,vtk_name,fig_size,bbox_to_anchor, local_range)
    cropping(path+data_name+'.png',data_name)
    # plot_3d_point(data_name, ttk_name, path,vtk_name)    
else:
    plot_fuc(data_name, ttk_name, path,ply_name,bbox_to_anchor,fig_size)

    # plot_3d_point(i)

