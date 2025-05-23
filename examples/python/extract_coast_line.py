from paraview.simple import *
import argparse


parser = argparse.ArgumentParser(description='rescale and extract one block from vtp.')
parser.add_argument('-p', '--path', type=str, default='path', help='path to all files')
parser.add_argument('-i', '--input', type=str, default='input.vtp', help='input vtp file')
parser.add_argument('-o', '--output', type=str, default='default.vtk', help='output vtp file')
parser.add_argument('-r', '--local_range', type=str, default='0,1,0,1,0,1', help='the range of points')


args = parser.parse_args()
# File paths
input_file = args.path + args.input
output_file = args.path + args.output
range =  tuple(map(float, args.local_range.split('-')))

# Load the VTP file
polyData = XMLPolyDataReader(FileName=input_file)

# Apply a transformation (Translate and then Scale)
transform = Transform(registrationName='Transform1',Input=polyData)
transform.Transform = 'Transform'

transform.Transform.Translate = [1800, 900, 0]  # Then translate
transform.Transform.Scale = [10, 10, 1]  # Scale the coordinates



# Update the pipeline to apply transformations
UpdatePipeline()


# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=transform)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['CELLS', 'GMI_CNTRY']
# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [(3600*range[0]+1800)%3600, 1800*range[2], 0.0]
clip1.HyperTreeGridClipper.Normal = [-1, 0, 0]
# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [(3600*range[0]+1800)%3600, 1800*range[2], 0.0]
clip1.ClipType.Normal = [-1, 0, 0]
# Update the pipeline to apply clipping
UpdatePipeline()

clip2 = Clip(registrationName='Clip2', Input=clip1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['CELLS', 'GMI_CNTRY']
# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [(3600*range[0]+1800)%3600, 1800*range[2], 0.0]
clip2.HyperTreeGridClipper.Normal = [0, -1, 0]
# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [(3600*range[0]+1800)%3600, 1800*range[2], 0.0]
clip2.ClipType.Normal = [0, -1, 0]
UpdatePipeline()

clip3 = Clip(registrationName='Clip2', Input=clip2)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['CELLS', 'GMI_CNTRY']
# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [(3600*range[1]+1800)%3600, 1800*range[3], 0.0]
clip3.HyperTreeGridClipper.Normal = [0, 1, 0]
# Properties modified on clip3.ClipType
clip3.ClipType.Origin = [(3600*range[1]+1800)%3600, 1800*range[3], 0.0]
clip3.ClipType.Normal = [0, 1, 0]
UpdatePipeline()

clip4 = Clip(registrationName='Clip2', Input=clip3)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = ['CELLS', 'GMI_CNTRY']
# init the 'Plane' selected for 'HyperTreeGridClipper'
clip4.HyperTreeGridClipper.Origin = [(3600*range[1]+1800)%3600, 1800*range[3], 0.0]
clip4.HyperTreeGridClipper.Normal = [1, 0, 0]
# Properties modified on clip4.ClipType
clip4.ClipType.Origin = [(3600*range[1]+1800)%3600, 1800*range[3], 0.0]
clip4.ClipType.Normal = [1, 0, 0]
UpdatePipeline()

transform2 = Transform(registrationName='Transform2',Input=clip4)
transform2.Transform = 'Transform'

transform2.Transform.Translate = [1800, 0, 0]  # Then translate
transform2.Transform.Scale = [1, 1, 1]  # Scale the coordinates
UpdatePipeline()
# Save the clipped data to a new VTP file
SaveData(output_file, proxy=transform2)

