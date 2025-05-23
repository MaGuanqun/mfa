from paraview.simple import *
import pandas as pd
import os
import csv
import argparse
import time

#use environment "mfa" for paraview+ttk

def compute_critical_point(path,name, csv_file):  
    
    _, extension = os.path.splitext(name)
    if extension.lower() == '.vtk':
        qmcpackvtk = LegacyVTKReader(registrationName='qmcpack.vtk', FileNames=[path+name])

        # create a new 'Tetrahedralize'
        calculator1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=qmcpackvtk)

    else:
        up_1ply = PLYReader(registrationName=name, FileNames=[path+name])
        calculator1 = Calculator(registrationName='Calculator1', Input=up_1ply)
        calculator1.Function = 'coordsZ'
        
    calculator1.UpdatePipeline()
    

    start_time = time.perf_counter()

    tTKScalarFieldCriticalPoints1 = TTKScalarFieldCriticalPoints(registrationName='TTKScalarFieldCriticalPoints1', Input=calculator1)
    

    SaveData("CriticalPoints.vtk", proxy=tTKScalarFieldCriticalPoints1)
    end_time =  time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Computing critical points needs {elapsed_time} seconds.")
    
    vtkData = OpenDataFile("CriticalPoints.vtk")
    # outputPort = servermanager.Fetch(vtkData)
    dataObject = servermanager.Fetch(vtkData)
    pointData = dataObject.GetPointData()
    # print(pointData)
    
    criticalTypeArray = pointData.GetArray("CriticalType")
    isOnBoundaryArray = pointData.GetArray("IsOnBoundary")    
    if extension.lower() == '.vtk':
        valueArray = pointData.GetArray("var0")
        
    positionArray = dataObject.GetPoints().GetData()
    
    # print(criticalTypeArray)
    
    numPoints = pointData.GetNumberOfTuples()

# Open a CSV file for writing
    if extension.lower() == '.vtk':
        with open(path+csv_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            # Write the header
            csvwriter.writerow(["CriticalType", "IsOnBoundary", "PositionX", "PositionY", "PositionZ", "var0"])
            # Write the data rows
            for i in range(numPoints):
                criticalType =int(criticalTypeArray.GetTuple(i)[0])
                # print(criticalType)
                isOnBoundary = int(isOnBoundaryArray.GetTuple(i)[0])
                var0 = valueArray.GetTuple(i)[0]
                position = positionArray.GetTuple(i)  # This retrieves the (x, y, z) tuple            
                # Write the current row to the CSV
                csvwriter.writerow([criticalType, isOnBoundary, position[0], position[1], position[2],var0])
    else:
        with open(path+csv_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            # Write the header
            csvwriter.writerow(["CriticalType", "IsOnBoundary", "PositionX", "PositionY", "PositionZ"])
            # Write the data rows
            for i in range(numPoints):
                criticalType =int(criticalTypeArray.GetTuple(i)[0])
                # print(criticalType)
                isOnBoundary = int(isOnBoundaryArray.GetTuple(i)[0])
                position = positionArray.GetTuple(i)  # This retrieves the (x, y, z) tuple            
                # Write the current row to the CSV
                csvwriter.writerow([criticalType, isOnBoundary, position[0], position[1], position[2]])

    
    if os.path.exists("CriticalPoints.vtk"):
        os.remove("CriticalPoints.vtk")
    
    

def pre_process_csv_file(csv_file):
    df = pd.read_csv(csv_file)
    # df = df.drop(columns=['ttkVertexScalarField','Result','Result_Order'])
    df = df[df['IsOnBoundary'] != 1]
    
    len1=len(df['CriticalType'])
    df = df[df['CriticalType'] < 4]
    len2=len(df['CriticalType'])
    
    # if(len1!=len2):
    #     print("size of critical point changed from ",len1," to ",len2) 
    # df['Points:2']=410
    df.to_csv(csv_file, index=False)
    print("size of critical point ",len(df['CriticalType']))



parser = argparse.ArgumentParser(description='TTK-critical points.')

parser.add_argument('-p', '--path', type=str, default='path', help='path to all files')
parser.add_argument('-i', '--input_name', type=str, default='file_name.ply', help='file to compute critical points')
parser.add_argument('-o','--output_name', type=str, default='ttk.csv', help='output csv name')


args = parser.parse_args()

path=args.path
file_name=args.input_name
csv_file=args.output_name
# path = '../../build/examples/scale/'
# file_name ='scale.ply'
# csv_file='ttk'+'.csv'
# file_name ='qmcpack.vtk'
# csv_file='ttk_'+'.csv'
# file_name='scale_'+str(i)+'_block_1.ply'
# csv_file='ttk_'+str(i)+'_block_1.csv'
compute_critical_point(path,file_name,csv_file )
pre_process_csv_file(path + csv_file)
print('Finish computing critical points for patch_'+file_name)
