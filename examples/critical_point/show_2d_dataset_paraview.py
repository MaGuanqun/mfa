# state file generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Light'
light1 = CreateLight()
light1.Intensity = 0.1
light1.Position = [1301.206157779645, 588.8039315980793, 1283.8938628886613]
light1.FocalPoint = [1259.9999999999998, 629.9999999999999, 1.0204630217819217e-13]

# create light
# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2950, 1568]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [1260.0, 630.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1301.206157779645, 588.8039315980793, 1283.8938628886613]
renderView1.CameraFocalPoint = [1259.9999999999998, 629.9999999999999, 1.0204630217819217e-13]
renderView1.CameraViewUp = [0.0012258075109972144, 0.9994861245922694, 0.032030987238494915]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 402.49223594996215
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.AdditionalLights = light1
renderView1.InteractionMode = '2D'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(2950, 1568)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------


num='1'

folder_name= 's3d'
model_name = folder_name
up=''



coast_line_file_name='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/coastLine'+num+'.vtk'
cpe_without_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/'+model_name+'_'+num+up+'_without_neighbor.csv'
cpe_with_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/'+model_name+'_'+num+up+'_with_neighbor.csv'
ttk_without_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/ttk_'+num+up+'_without_neighbor.csv'
block_ply_file='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/'+model_name+'_'+num+up+'.ply'
ttk_with_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/ttk_'+num+up+'_with_neighbor.csv'
ttk_result='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/ttk_'+num+up+'.csv'
cesm_result='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+folder_name+'/'+model_name+'_'+num+up+'.csv'


# create a new 'Legacy VTK Reader'
if model_name == 'cesm':
    coastLine2vtk = LegacyVTKReader(registrationName='coastLine2.vtk', FileNames=[coast_line_file_name])
    color_map_max = 430
    radius=4
if model_name == 's3d':
    color_map_max = 360
    if int(num)>3:
        radius = 1.0
    else:
        radius=1.5

# create a new 'CSV Reader'
cesm_2_without_neighborcsv = CSVReader(registrationName='cesm_2_without_neighbor.csv', FileName=[cpe_without_neighbor])


# create a new 'CSV Reader'
cesm_2_with_neighborcsv = CSVReader(registrationName='cesm_2_with_neighbor.csv', FileName=[cpe_with_neighbor])

# create a new 'CSV Reader'
ttk_2_without_neighborcsv = CSVReader(registrationName='ttk_2_without_neighbor.csv', FileName=[ttk_without_neighbor])

# create a new 'PLY Reader'
cesm_2ply = PLYReader(registrationName='cesm_2.ply', FileNames=[block_ply_file])

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=cesm_2ply)
calculator1.ResultArrayName = 'oriZ'
calculator1.Function = 'coordsZ'

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
calculator2.CoordinateResults = 1
calculator2.ResultArrayName = 'FlatZ'
calculator2.Function = 'coordsX*iHat + coordsY*jHat + 0*kHat'

# create a new 'Table To Points'
tableToPoints4 = TableToPoints(registrationName='TableToPoints4', Input=cesm_2_with_neighborcsv)
tableToPoints4.XColumn = 'x0'
tableToPoints4.YColumn = 'x1'
tableToPoints4.ZColumn = 'x0'
tableToPoints4.a2DPoints = 1

# create a new 'CSV Reader'
cesm_2csv = CSVReader(registrationName='cesm_2.csv', FileName=[cesm_result])

# # create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=cesm_2csv)
tableToPoints1.XColumn = 'x0'
tableToPoints1.YColumn = 'x1'
tableToPoints1.ZColumn = 'x2'
tableToPoints1.a2DPoints = 1

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints1 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints1', Input=tableToPoints1)
tTKIcospheresFromPoints1.Radius = radius


# show data from tTKIcospheresFromPoints3
tTKIcospheresFromPoints1Display = Show(tTKIcospheresFromPoints1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints1Display.Representation = 'Surface'
tTKIcospheresFromPoints1Display.AmbientColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints1Display.ColorArrayName = [None, '']
tTKIcospheresFromPoints1Display.DiffuseColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints1Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints1Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints1Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints1Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints1Display.ScaleFactor = 70.425
tTKIcospheresFromPoints1Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints1Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints1Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints1Display.GaussianRadius = 3.52125
tTKIcospheresFromPoints1Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints1Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints1Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints1Display.WriteLog = ''
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]


# create a new 'Table To Points'
tableToPoints3 = TableToPoints(registrationName='TableToPoints3', Input=cesm_2_without_neighborcsv)
tableToPoints3.XColumn = 'x0'
tableToPoints3.YColumn = 'x1'
tableToPoints3.ZColumn = 'x0'
tableToPoints3.a2DPoints = 1

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints4 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints4', Input=tableToPoints3)
tTKIcospheresFromPoints4.Radius = radius

# create a new 'PLY Reader'
# cesm_1ply = PLYReader(registrationName='cesm_1.ply', FileNames=['//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/cesm/cesm_1.ply'])

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints3 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints3', Input=tableToPoints4)
tTKIcospheresFromPoints3.Radius = radius

# create a new 'CSV Reader'
ttk_2_with_neighborcsv = CSVReader(registrationName='ttk_2_with_neighbor.csv', FileName=[ttk_with_neighbor])

# create a new 'Table To Points'
tableToPoints5 = TableToPoints(registrationName='TableToPoints5', Input=ttk_2_with_neighborcsv)
tableToPoints5.XColumn = 'x0'
tableToPoints5.YColumn = 'x1'
tableToPoints5.ZColumn = 'x0'
tableToPoints5.a2DPoints = 1

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints5 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints5', Input=tableToPoints5)
tTKIcospheresFromPoints5.Radius = radius

# create a new 'Table To Points'
tableToPoints6 = TableToPoints(registrationName='TableToPoints6', Input=ttk_2_without_neighborcsv)
tableToPoints6.XColumn = 'x0'
tableToPoints6.YColumn = 'x1'
tableToPoints6.ZColumn = 'x0'
tableToPoints6.a2DPoints = 1

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints6 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints6', Input=tableToPoints6)
tTKIcospheresFromPoints6.Radius = radius

# create a new 'CSV Reader'
ttk_2csv = CSVReader(registrationName='ttk_2.csv', FileName=[ttk_result])

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(registrationName='TableToPoints2', Input=ttk_2csv)
tableToPoints2.XColumn = 'PositionX'
tableToPoints2.YColumn = 'PositionY'
tableToPoints2.ZColumn = 'CriticalType'
tableToPoints2.a2DPoints = 1

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints2 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints2', Input=tableToPoints2)
tTKIcospheresFromPoints2.Radius = radius


# show data from tTKIcospheresFromPoints6
tTKIcospheresFromPoints2Display = Show(tTKIcospheresFromPoints2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints2Display.Representation = 'Surface'
tTKIcospheresFromPoints2Display.AmbientColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints2Display.ColorArrayName = [None, '']
tTKIcospheresFromPoints2Display.DiffuseColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints2Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints2Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints2Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints2Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints2Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints2Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints2Display.ScaleFactor = 65.2
tTKIcospheresFromPoints2Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints2Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints2Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints2Display.GaussianRadius = 3.2600000000000002
tTKIcospheresFromPoints2Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints2Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints2Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints2Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints2Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints2Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints2Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints2Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints2Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]



# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

if model_name == 'cesm':
    # show data from coastLine2vtk
    coastLine2vtkDisplay = Show(coastLine2vtk, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    coastLine2vtkDisplay.Representation = 'Wireframe'
    coastLine2vtkDisplay.AmbientColor = [0.5254901960784314, 0.5254901960784314, 0.5254901960784314]
    coastLine2vtkDisplay.ColorArrayName = [None, '']
    coastLine2vtkDisplay.DiffuseColor = [0.5254901960784314, 0.5254901960784314, 0.5254901960784314]
    coastLine2vtkDisplay.LineWidth = 4.0
    coastLine2vtkDisplay.SelectTCoordArray = 'None'
    coastLine2vtkDisplay.SelectNormalArray = 'None'
    coastLine2vtkDisplay.SelectTangentArray = 'None'
    coastLine2vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    coastLine2vtkDisplay.SelectOrientationVectors = 'None'
    coastLine2vtkDisplay.ScaleFactor = 49.09169921875
    coastLine2vtkDisplay.SelectScaleArray = 'None'
    coastLine2vtkDisplay.GlyphType = 'Arrow'
    coastLine2vtkDisplay.GlyphTableIndexArray = 'None'
    coastLine2vtkDisplay.GaussianRadius = 2.4545849609375
    coastLine2vtkDisplay.SetScaleArray = [None, '']
    coastLine2vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    coastLine2vtkDisplay.OpacityArray = [None, '']
    coastLine2vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    coastLine2vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
    coastLine2vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
    coastLine2vtkDisplay.ScalarOpacityUnitDistance = 44.81089085482233
    coastLine2vtkDisplay.OpacityArrayName = ['CELLS', 'GMI_CNTRY']
    coastLine2vtkDisplay.SelectInputVectors = [None, '']
    coastLine2vtkDisplay.WriteLog = ''

# show data from calculator2
calculator2Display = Show(calculator2, renderView1, 'GeometryRepresentation')
calculator2Display.SetScalarBarVisibility(renderView1, False)
# get 2D transfer function for 'oriZ'
oriZTF2D = GetTransferFunction2D('oriZ')
oriZTF2D.ScalarRangeInitialized = 1
oriZTF2D.Range = [-1.0, color_map_max, 4.684686660766602, 419.4297790527344]

# get color transfer function/color map for 'oriZ'
oriZLUT = GetColorTransferFunction('oriZ')
oriZLUT.AutomaticRescaleRangeMode = 'Never'
oriZLUT.TransferFunction2D = oriZTF2D
oriZLUT.RGBPoints = [-1.0, 0.231373, 0.298039, 0.752941, 107.34529594091914, 0.5529411764705883, 0.6901960784313725, 0.996078431372549, 214.5, 0.865003, 0.865003, 0.865003, color_map_max, 0.705882, 0.0156863, 0.14902]
oriZLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.

calculator2Display.Representation = 'Surface'
calculator2Display.ColorArrayName = ['POINTS', 'oriZ']
calculator2Display.LookupTable = oriZLUT
calculator2Display.SelectTCoordArray = 'None'
calculator2Display.SelectNormalArray = 'None'
calculator2Display.SelectTangentArray = 'None'
calculator2Display.OSPRayScaleArray = 'oriZ'
calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator2Display.SelectOrientationVectors = 'FlatZ'
calculator2Display.ScaleFactor = 72.0
calculator2Display.SelectScaleArray = 'oriZ'
calculator2Display.GlyphType = 'Arrow'
calculator2Display.GlyphTableIndexArray = 'oriZ'
calculator2Display.GaussianRadius = 3.6
calculator2Display.SetScaleArray = ['POINTS', 'oriZ']
calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display.OpacityArray = ['POINTS', 'oriZ']
calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
calculator2Display.PolarAxes = 'PolarAxesRepresentation'
calculator2Display.SelectInputVectors = ['POINTS', 'FlatZ']
calculator2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator2Display.ScaleTransferFunction.Points = [268.31219482421875, 0.0, 0.5, 0.0, 419.4297790527344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator2Display.OpacityTransferFunction.Points = [268.31219482421875, 0.0, 0.5, 0.0, 419.4297790527344, 1.0, 0.5, 0.0]

# show data from tTKIcospheresFromPoints3
tTKIcospheresFromPoints3Display = Show(tTKIcospheresFromPoints3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints3Display.Representation = 'Surface'
tTKIcospheresFromPoints3Display.AmbientColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints3Display.ColorArrayName = [None, '']
tTKIcospheresFromPoints3Display.DiffuseColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints3Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints3Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints3Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints3Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints3Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints3Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints3Display.ScaleFactor = 70.425
tTKIcospheresFromPoints3Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints3Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints3Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints3Display.GaussianRadius = 3.52125
tTKIcospheresFromPoints3Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints3Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints3Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints3Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints3Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints3Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints3Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints3Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints3Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints3Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from tTKIcospheresFromPoints4
tTKIcospheresFromPoints4Display = Show(tTKIcospheresFromPoints4, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints4Display.Representation = 'Surface'

# change solid color
tTKIcospheresFromPoints4Display.AmbientColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints4Display.ColorArrayName = [None, '']
tTKIcospheresFromPoints4Display.DiffuseColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints4Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints4Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints4Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints4Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints4Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints4Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints4Display.ScaleFactor = 72.16890000000001
tTKIcospheresFromPoints4Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints4Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints4Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints4Display.GaussianRadius = 3.6084450000000006
tTKIcospheresFromPoints4Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints4Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints4Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints4Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints4Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints4Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints4Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints4Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints4Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints4Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from tTKIcospheresFromPoints5
tTKIcospheresFromPoints5Display = Show(tTKIcospheresFromPoints5, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints5Display.Representation = 'Surface'
tTKIcospheresFromPoints5Display.AmbientColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints5Display.ColorArrayName = [None, '']
tTKIcospheresFromPoints5Display.DiffuseColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints5Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints5Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints5Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints5Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints5Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints5Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints5Display.ScaleFactor = 70.4
tTKIcospheresFromPoints5Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints5Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints5Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints5Display.GaussianRadius = 3.52
tTKIcospheresFromPoints5Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints5Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints5Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints5Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints5Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints5Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints5Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints5Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints5Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints5Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from tTKIcospheresFromPoints6
tTKIcospheresFromPoints6Display = Show(tTKIcospheresFromPoints6, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints6Display.Representation = 'Surface'
tTKIcospheresFromPoints6Display.AmbientColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints6Display.ColorArrayName = [None, '']
tTKIcospheresFromPoints6Display.DiffuseColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints6Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints6Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints6Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints6Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints6Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints6Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints6Display.ScaleFactor = 65.2
tTKIcospheresFromPoints6Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints6Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints6Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints6Display.GaussianRadius = 3.2600000000000002
tTKIcospheresFromPoints6Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints6Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints6Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints6Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints6Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints6Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints6Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints6Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints6Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints6Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for oriZLUT in view renderView1
oriZLUTColorBar = GetScalarBar(oriZLUT, renderView1)
oriZLUTColorBar.WindowLocation = 'Any Location'
oriZLUTColorBar.Position = [0.8624152542372882, 0.3571428571428572]
oriZLUTColorBar.Title = 'oriZ'
oriZLUTColorBar.ComponentTitle = ''
oriZLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
oriZLUTColorBar.ScalarBarLength = 0.33000000000000007

# set color bar visibility
oriZLUTColorBar.Visibility = 0

# show color legend
calculator2Display.SetScalarBarVisibility(renderView1, False)


renderView1 = GetActiveViewOrCreate('RenderView')
Hide(tTKIcospheresFromPoints2, renderView1)
Hide(tTKIcospheresFromPoints1, renderView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'oriZ'
oriZPWF = GetOpacityTransferFunction('oriZ')
oriZPWF.Points = [-1.0, 0.0, 0.5, 0.0, color_map_max, 1.0, 0.5, 0.0]
oriZPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(tTKIcospheresFromPoints4)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')