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



model_name='rti'
num='3'
up=''

block_num=4

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2950, 1582]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraFocalDisk = 1.0
if model_name == 'rti':
    if num == '2':
        # Create a new 'Render View'
        renderView1.CenterOfRotation = [64.92785, 29.0747, 160.98700000000002]
        renderView1.CameraPosition = [63.38366693949249, 9.59415223521492, 200.14904080222675]
        renderView1.CameraFocalPoint = [64.92785000000026, 29.07469999999993, 160.98700000000034]
        renderView1.CameraViewUp = [0.999361485268149, -0.010661062872572851, 0.034102250674497646]
        renderView1.CameraParallelScale = 11.327707698934507
    if num == '1':
        renderView1.CenterOfRotation = [81.9973, 72.00165000000001, 93.01920000000001]
        renderView1.CameraPosition = [89.34871163410003, 110.3832835382401, 111.52986455133653]
        renderView1.CameraFocalPoint = [81.99729999999998, 72.00164999999964, 93.01920000000025]
        renderView1.CameraViewUp = [0.9854245994502193, -0.15049225920260215, -0.07931228604996649]
        renderView1.CameraParallelScale = 11.191754940691839
    if num =='3':
        renderView1.CenterOfRotation = [52.07415, 56.99205, 209.0125]
        renderView1.CameraPosition = [54.286142022594476, 96.16470452607393, 228.5821035297418]
        renderView1.CameraFocalPoint = [52.07415, 56.99205, 209.0125]
        renderView1.CameraViewUp = [0.9987259273067689, -0.045634273432460265, -0.0215414765864826]
        renderView1.CameraParallelScale = 11.347847623888862
if model_name == 'qmcpack':
    if num == '1':
        renderView1.CameraPosition = [4.576380357038572, 123.92906419556861, 52.21384640436362]
        renderView1.CameraFocalPoint = [58.058449999999944, 58.01374999999994, 63.000000000000064]
        renderView1.CameraViewUp = [-0.0763681799562973, 0.10036668355129398, 0.9920153375442726]
        renderView1.CenterOfRotation = [58.05845, 58.01375, 63.0]
        renderView1.CameraParallelScale = 22.14605169019977
    if num =='2':
        renderView1.CameraPosition = [16.08254620450559, -50.839407456730584, 79.28367655247452]
        renderView1.CameraFocalPoint = [34.0, 34.0, 85.5]
        renderView1.CameraViewUp = [-0.012324623933711224, -0.07048098521759592, 0.9974369826548692]
        renderView1.CenterOfRotation = [34.0, 34.0, 85.5]
    if  num == '3':
        renderView1.CameraPosition = [92.24402703013402, -52.26137788076484, 79.49065345253567]
        renderView1.CameraFocalPoint = [50.99999999999995, 17.02520000000006, 51.00000000000003]
        renderView1.CameraViewUp = [-0.16474327362572205, 0.28963494385290084, 0.9428527207864479]
        renderView1.CameraParallelScale = 22.133821248035776
        renderView1.CenterOfRotation = [51.0, 17.025199999999998, 51.0]


renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
SetActiveView(None)


# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(2950, 1582)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

rti_ply_file='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+num+'.ply'
ttk_without_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/ttk_'+num+up+'_without_neighbor.csv'
ttk_result='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/ttk_'+num+up+'.csv'
rti_without_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+num+up+'_without_neighbor.csv'
rti_result='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+num+up+'.csv'
ttk_with_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/ttk_'+num+up+'_with_neighbor.csv'
rti_with_neighbor='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+num+up+'_with_neighbor.csv'
rti_model='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+num+up+'.vtk'

rti_ori_model='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+num+'.vtk'
block_ply_file='//wsl.localhost/Ubuntu/home/guanqunma/mfa/build/examples/'+model_name+'/'+model_name+'_'+str(block_num)+'.ply'


if model_name == 'rti':
    radius = 0.15
    if num == '2' or num =='5':
        lower_threshold=8
        higher_threshold=16
        iso_surface_value=13.5
    if num == '1' or num =='4':
        lower_threshold = 2
        higher_threshold= 5
        iso_surface_value = 3.5
    if num == '3' or num =='6':
        lower_threshold = 8
        higher_threshold = 13
        iso_surface_value = 13
if model_name == 'qmcpack':
    radius = 0.4
    if num == '1' or num =='4':
        lower_threshold=-0.000006352
        higher_threshold=0.001
        iso_surface_value=0.00002
    if num == '2' or num == '5':
        lower_threshold=-0.000024488
        higher_threshold=0.002
        iso_surface_value=0.0005  
    if num == '3' or num =='6':
        lower_threshold=-6.29358e-06
        higher_threshold=6.29358e-05
        iso_surface_value=-6.29358e-06

if up=='_up':
    if num=='3' or num =='1' or num=='2':
        qmcpack_4ply = PLYReader(registrationName='qmcpack_4.ply', FileNames=[block_ply_file])

# create a new 'PLY Reader'
rti_2ply = PLYReader(registrationName='rti_2.ply', FileNames=[rti_ply_file])

# create a new 'CSV Reader'
ttk_2_up_without_neighborcsv = CSVReader(registrationName='ttk_2_up_without_neighbor.csv', FileName=[ttk_without_neighbor])

# create a new 'CSV Reader'
ttk_2_upcsv = CSVReader(registrationName='ttk_2_up.csv', FileName=[ttk_result])

# create a new 'CSV Reader'
rti_2_up_without_neighborcsv = CSVReader(registrationName='rti_2_up_without_neighbor.csv', FileName=[rti_without_neighbor])

# create a new 'Table To Points'
tableToPoints4 = TableToPoints(registrationName='TableToPoints4', Input=rti_2_up_without_neighborcsv)
tableToPoints4.XColumn = 'x0'
tableToPoints4.YColumn = 'x1'
tableToPoints4.ZColumn = 'x2'

# create a new 'CSV Reader'
rti_2_upcsv = CSVReader(registrationName='rti_2_up.csv', FileName=[rti_result])


# create a new 'Table To Points'
tableToPoints2 = TableToPoints(registrationName='TableToPoints2', Input=rti_2_upcsv)
tableToPoints2.XColumn = 'x0'
tableToPoints2.YColumn = 'x1'
tableToPoints2.ZColumn = 'x2'

# create a new 'Tetrahedralize'
tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=tableToPoints2)

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints8 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints8', Input=tableToPoints2)
tTKIcospheresFromPoints8.Radius = 0.075

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=tetrahedralize1)
threshold1.Scalars = ['POINTS', 'x3']
threshold1.LowerThreshold = lower_threshold
threshold1.UpperThreshold = higher_threshold

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints2 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints2', Input=threshold1)
tTKIcospheresFromPoints2.Radius = radius



# show data from rti_2ply
tTKIcospheresFromPoints2Display = Show(tTKIcospheresFromPoints2, renderView1, 'GeometryRepresentation')

tTKIcospheresFromPoints2Display.Representation = 'Surface'
tTKIcospheresFromPoints2Display.AmbientColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints2Display.ColorArrayName = ['POINTS', '']
tTKIcospheresFromPoints2Display.DiffuseColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints2Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints2Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints2Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints2Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints2Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints2Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints2Display.ScaleFactor = 1.5650003814697273
tTKIcospheresFromPoints2Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints2Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints2Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints2Display.GaussianRadius = 0.07825001907348636
tTKIcospheresFromPoints2Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints2Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints2Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints2Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints2Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints2Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints2Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints2Display.WriteLog = ''

Hide(tTKIcospheresFromPoints2, renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=ttk_2_upcsv)
tableToPoints1.XColumn = 'PositionX'
tableToPoints1.YColumn = 'PositionY'
tableToPoints1.ZColumn = 'PositionZ'

# create a new 'Tetrahedralize'
tetrahedralize5 = Tetrahedralize(registrationName='Tetrahedralize5', Input=tableToPoints4)

# create a new 'CSV Reader'
ttk_2_up_with_neighborcsv = CSVReader(registrationName='ttk_2_up_with_neighbor.csv', FileName=[ttk_with_neighbor])

# create a new 'Table To Points'
tableToPoints5 = TableToPoints(registrationName='TableToPoints5', Input=ttk_2_up_with_neighborcsv)
tableToPoints5.XColumn = 'x0'
tableToPoints5.YColumn = 'x1'
tableToPoints5.ZColumn = 'x2'

# create a new 'Tetrahedralize'
tetrahedralize4 = Tetrahedralize(registrationName='Tetrahedralize4', Input=tableToPoints5)

# create a new 'Table To Points'
tableToPoints6 = TableToPoints(registrationName='TableToPoints6', Input=ttk_2_up_without_neighborcsv)
tableToPoints6.XColumn = 'x0'
tableToPoints6.YColumn = 'x1'
tableToPoints6.ZColumn = 'x2'

# create a new 'Tetrahedralize'
tetrahedralize3 = Tetrahedralize(registrationName='Tetrahedralize3', Input=tableToPoints6)

# create a new 'Threshold'
threshold6 = Threshold(registrationName='Threshold6', Input=tetrahedralize3)
threshold6.Scalars = ['POINTS', 'x3']
threshold6.LowerThreshold = lower_threshold
threshold6.UpperThreshold = higher_threshold

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints3 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints3', Input=threshold6)
tTKIcospheresFromPoints3.Radius = radius

# create a new 'Threshold'
threshold4 = Threshold(registrationName='Threshold4', Input=tetrahedralize5)
threshold4.Scalars = ['POINTS', 'x3']
threshold4.LowerThreshold = lower_threshold
threshold4.UpperThreshold = higher_threshold

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints5 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints5', Input=threshold4)
tTKIcospheresFromPoints5.Radius = radius

# create a new 'CSV Reader'
rti_2_up_with_neighborcsv = CSVReader(registrationName='rti_2_up_with_neighbor.csv', FileName=[rti_with_neighbor])

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(registrationName='TableToPoints3', Input=rti_2_up_with_neighborcsv)
tableToPoints3.XColumn = 'x0'
tableToPoints3.YColumn = 'x1'
tableToPoints3.ZColumn = 'x2'

# create a new 'Tetrahedralize'
tetrahedralize6 = Tetrahedralize(registrationName='Tetrahedralize6', Input=tableToPoints3)

# create a new 'Threshold'
threshold3 = Threshold(registrationName='Threshold3', Input=tetrahedralize6)
threshold3.Scalars = ['POINTS', 'x3']
threshold3.LowerThreshold = lower_threshold
threshold3.UpperThreshold = higher_threshold

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints6 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints6', Input=threshold3)
tTKIcospheresFromPoints6.Radius = radius

# create a new 'Threshold'
threshold5 = Threshold(registrationName='Threshold5', Input=tetrahedralize4)
threshold5.Scalars = ['POINTS', 'x3']
threshold5.LowerThreshold = lower_threshold
threshold5.UpperThreshold = higher_threshold

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints4 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints4', Input=threshold5)
tTKIcospheresFromPoints4.Radius = radius


if int(num) < 4 or model_name =='rti':
    rti_2_upvtk = LegacyVTKReader(registrationName='rti_2_up.vtk', FileNames=[rti_model])
    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=rti_2_upvtk)
    contour1.ContourBy = ['POINTS', 'var0']
    contour1.Isosurfaces = [iso_surface_value]
    contour1.PointMergeMethod = 'Uniform Binning'

if int(num)<4 and model_name == 'rti':
    isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=rti_2_upvtk)
    isoVolume1.InputScalars = ['POINTS', 'var0']
    isoVolume1.ThresholdRange = [lower_threshold, higher_threshold]
elif up=='_up':
    rti_ori_upvtk = LegacyVTKReader(registrationName='rti_2_ori.vtk', FileNames=[rti_ori_model])
    # create a new 'Iso Volume'
    isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=rti_ori_upvtk)
    isoVolume1.InputScalars = ['POINTS', 'var0']
    isoVolume1.ThresholdRange = [lower_threshold, higher_threshold]
elif int(num) < 4:
    # create a new 'Legacy VTK Reader'
    # create a new 'Iso Volume'
    isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=rti_2_upvtk)
    isoVolume1.InputScalars = ['POINTS', 'var0']
    isoVolume1.ThresholdRange = [lower_threshold, higher_threshold]



# create a new 'Tetrahedralize'
tetrahedralize2 = Tetrahedralize(registrationName='Tetrahedralize2', Input=tableToPoints1)

# create a new 'Threshold'
threshold2 = Threshold(registrationName='Threshold2', Input=tetrahedralize2)
threshold2.Scalars = ['POINTS', 'var0']
threshold2.LowerThreshold = lower_threshold
threshold2.UpperThreshold = higher_threshold

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints1 = TTKIcospheresFromPoints(registrationName='TTKIcospheresFromPoints1', Input=threshold2)
tTKIcospheresFromPoints1.Radius = radius


# show data from rti_2ply
tTKIcospheresFromPoints1Display = Show(tTKIcospheresFromPoints1, renderView1, 'GeometryRepresentation')

tTKIcospheresFromPoints1Display.Representation = 'Surface'
tTKIcospheresFromPoints1Display.AmbientColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints1Display.ColorArrayName = ['POINTS', '']
tTKIcospheresFromPoints1Display.DiffuseColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints1Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints1Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints1Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints1Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints1Display.ScaleFactor = 1.5650003814697273
tTKIcospheresFromPoints1Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints1Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints1Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints1Display.GaussianRadius = 0.07825001907348636
tTKIcospheresFromPoints1Display.SetScaleArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.OpacityArray = ['POINTS', 'Normals']
tTKIcospheresFromPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints1Display.PolarAxes = 'PolarAxesRepresentation'
tTKIcospheresFromPoints1Display.SelectInputVectors = ['POINTS', 'Normals']
tTKIcospheresFromPoints1Display.WriteLog = ''


Hide(tTKIcospheresFromPoints1, renderView1)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from rti_2ply
rti_2plyDisplay = Show(rti_2ply, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
rti_2plyDisplay.Representation = 'Wireframe'
rti_2plyDisplay.AmbientColor = [0.6392156862745098, 0.6392156862745098, 0.6392156862745098]
rti_2plyDisplay.ColorArrayName = ['POINTS', '']
rti_2plyDisplay.DiffuseColor = [0.6392156862745098, 0.6392156862745098, 0.6392156862745098]
rti_2plyDisplay.LineWidth = 3.0
rti_2plyDisplay.SelectTCoordArray = 'None'
rti_2plyDisplay.SelectNormalArray = 'Normals'
rti_2plyDisplay.SelectTangentArray = 'None'
rti_2plyDisplay.OSPRayScaleArray = 'Normals'
rti_2plyDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
rti_2plyDisplay.SelectOrientationVectors = 'None'
rti_2plyDisplay.ScaleFactor = 1.6
rti_2plyDisplay.SelectScaleArray = 'None'
rti_2plyDisplay.GlyphType = 'Arrow'
rti_2plyDisplay.GlyphTableIndexArray = 'None'
rti_2plyDisplay.GaussianRadius = 0.08
rti_2plyDisplay.SetScaleArray = ['POINTS', 'Normals']
rti_2plyDisplay.ScaleTransferFunction = 'PiecewiseFunction'
rti_2plyDisplay.OpacityArray = ['POINTS', 'Normals']
rti_2plyDisplay.OpacityTransferFunction = 'PiecewiseFunction'
rti_2plyDisplay.DataAxesGrid = 'GridAxesRepresentation'
rti_2plyDisplay.PolarAxes = 'PolarAxesRepresentation'
rti_2plyDisplay.SelectInputVectors = ['POINTS', 'Normals']
rti_2plyDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
rti_2plyDisplay.ScaleTransferFunction.Points = [-0.5773502588272095, 0.0, 0.5, 0.0, 0.5773502588272095, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
rti_2plyDisplay.OpacityTransferFunction.Points = [-0.5773502588272095, 0.0, 0.5, 0.0, 0.5773502588272095, 1.0, 0.5, 0.0]

# show data from tTKIcospheresFromPoints3
tTKIcospheresFromPoints3Display = Show(tTKIcospheresFromPoints3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tTKIcospheresFromPoints3Display.Representation = 'Surface'
tTKIcospheresFromPoints3Display.AmbientColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints3Display.ColorArrayName = ['POINTS', '']
tTKIcospheresFromPoints3Display.DiffuseColor = [0.4, 0.27450980392156865, 0.6196078431372549]
tTKIcospheresFromPoints3Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints3Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints3Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints3Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints3Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints3Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints3Display.ScaleFactor = 1.5650003814697273
tTKIcospheresFromPoints3Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints3Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints3Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints3Display.GaussianRadius = 0.07825001907348636
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
tTKIcospheresFromPoints4Display.AmbientColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints4Display.ColorArrayName = ['POINTS', '']
tTKIcospheresFromPoints4Display.DiffuseColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints4Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints4Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints4Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints4Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints4Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints4Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints4Display.ScaleFactor = 1.5749998474121107
tTKIcospheresFromPoints4Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints4Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints4Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints4Display.GaussianRadius = 0.07874999237060554
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
tTKIcospheresFromPoints5Display.AmbientColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints5Display.ColorArrayName = ['POINTS', '']
tTKIcospheresFromPoints5Display.DiffuseColor = [1.0, 0.3254901960784314, 0.8745098039215686]
tTKIcospheresFromPoints5Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints5Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints5Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints5Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints5Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints5Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints5Display.ScaleFactor = 1.1848400000000006
tTKIcospheresFromPoints5Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints5Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints5Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints5Display.GaussianRadius = 0.059242000000000024
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
tTKIcospheresFromPoints6Display.AmbientColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints6Display.ColorArrayName = ['POINTS', '']
tTKIcospheresFromPoints6Display.DiffuseColor = [1.0, 1.0, 0.0]
tTKIcospheresFromPoints6Display.SelectTCoordArray = 'None'
tTKIcospheresFromPoints6Display.SelectNormalArray = 'Normals'
tTKIcospheresFromPoints6Display.SelectTangentArray = 'None'
tTKIcospheresFromPoints6Display.OSPRayScaleArray = 'Normals'
tTKIcospheresFromPoints6Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints6Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints6Display.ScaleFactor = 1.58264
tTKIcospheresFromPoints6Display.SelectScaleArray = 'None'
tTKIcospheresFromPoints6Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints6Display.GlyphTableIndexArray = 'None'
tTKIcospheresFromPoints6Display.GaussianRadius = 0.079132
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


    # show data from isoVolume1
isoVolume1Display = Show(isoVolume1, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'var0'
var0TF2D = GetTransferFunction2D('var0')
var0TF2D.ScalarRangeInitialized = 1
var0TF2D.Range = [-10.0, 35.0, 4.817086219787598, 40]

# get color transfer function/color map for 'var0'
var0LUT = GetColorTransferFunction('var0')
var0LUT.AutomaticRescaleRangeMode = 'Never'
var0LUT.TransferFunction2D = var0TF2D
var0LUT.RGBPoints = [-10, 0.47, 0.84, 1.0, 4.685082912445069, 0.47, 0.84, 1.0, 19.5, 0.47, 0.84, 1.0, 40.0, 0.47, 0.84, 1.0]
var0LUT.ScalarRangeInitialized = 1.0


# get opacity transfer function/opacity map for 'var0'
var0PWF = GetOpacityTransferFunction('var0')

if model_name=='rti' and (num=='1' or num == '4'):
    var0PWF.Points = [4.0, 0.0, 0.5, 0.0, 35.0, 0.2, 0.5, 0.0]
elif  model_name=='rti' and (num=='2' or num == '5'):
    var0PWF.Points = [4.0, 0.0, 0.5, 0.0, 35.0, 0.3, 0.5, 0.0]
else:
    var0PWF.Points = [4.0, 0.0, 0.5, 0.0, 35.0, 0.5, 0.5, 0.0]
var0PWF.ScalarRangeInitialized = 1




# rescale color and/or opacity maps used to exactly fit the current data range
isoVolume1Display.RescaleTransferFunctionToDataRange(False, True)

# get color transfer function/color map for 'var0'
var0LUT = GetColorTransferFunction('var0')

# get opacity transfer function/opacity map for 'var0'
var0PWF = GetOpacityTransferFunction('var0')

# get 2D transfer function for 'var0'
var0TF2D = GetTransferFunction2D('var0')


if model_name=='qmcpack':
    if num == '2':
        var0PWF.Points = [-0.00012448809866327792, 0.0, 0.5, 0.0, 0.00014000007649883628, 0.05, 0.5, 0.0]
    if num =='3':
        var0PWF.Points = [-0.00012448809866327792, 0.0, 0.5, 0.0, 0.00014000007649883628, 0.04, 0.5, 0.0]


# trace defaults for the display properties.
isoVolume1Display.Representation = 'Volume'
isoVolume1Display.ColorArrayName = ['POINTS', 'var0']
isoVolume1Display.LookupTable = var0LUT
isoVolume1Display.SelectTCoordArray = 'None'
isoVolume1Display.SelectNormalArray = 'None'
isoVolume1Display.SelectTangentArray = 'None'
isoVolume1Display.OSPRayScaleArray = 'var0'
isoVolume1Display.OSPRayScaleFunction = 'PiecewiseFunction'
isoVolume1Display.SelectOrientationVectors = 'None'
isoVolume1Display.ScaleFactor = 1.6
isoVolume1Display.SelectScaleArray = 'var0'
isoVolume1Display.GlyphType = 'Arrow'
isoVolume1Display.GlyphTableIndexArray = 'var0'
isoVolume1Display.GaussianRadius = 0.08
isoVolume1Display.SetScaleArray = ['POINTS', 'var0']
isoVolume1Display.ScaleTransferFunction = 'PiecewiseFunction'
isoVolume1Display.OpacityArray = ['POINTS', 'var0']
isoVolume1Display.OpacityTransferFunction = 'PiecewiseFunction'
isoVolume1Display.DataAxesGrid = 'GridAxesRepresentation'
isoVolume1Display.PolarAxes = 'PolarAxesRepresentation'
isoVolume1Display.ScalarOpacityFunction = var0PWF
isoVolume1Display.ScalarOpacityUnitDistance = 0.2566556609830306
isoVolume1Display.OpacityArrayName = ['POINTS', 'var0']
isoVolume1Display.SelectInputVectors = ['POINTS', '']
isoVolume1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isoVolume1Display.ScaleTransferFunction.Points = [10.0, 0.0, 0.5, 0.0, 23.884660720825195, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isoVolume1Display.OpacityTransferFunction.Points = [10.0, 0.0, 0.5, 0.0, 23.884660720825195, 1.0, 0.5, 0.0]
isoVolume1Display.SetScalarBarVisibility(renderView1, False)

if int(num) < 4 or model_name =='rti':
    # show data from contour1
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.AmbientColor = [0.47058823529411764, 0.8431372549019608, 1.0]
    contour1Display.ColorArrayName = ['POINTS', '']
    contour1Display.DiffuseColor = [0.47058823529411764, 0.8431372549019608, 1.0]
    contour1Display.LookupTable = var0LUT
    contour1Display.SelectTCoordArray = 'None'
    contour1Display.SelectNormalArray = 'Normals'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = 'var0'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 1.6
    contour1Display.SelectScaleArray = 'var0'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'var0'
    contour1Display.GaussianRadius = 0.08
    contour1Display.SetScaleArray = ['POINTS', 'var0']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'var0']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    contour1Display.SelectInputVectors = ['POINTS', 'Normals']
    contour1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [12.5, 0.0, 0.5, 0.0, 12.501953125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [12.5, 0.0, 0.5, 0.0, 12.501953125, 1.0, 0.5, 0.0]

if up=='_up':
    if num =='1' or num == '2' or num =='3':
        # show data from qmcpack_4ply
        qmcpack_4plyDisplay = Show(qmcpack_4ply, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        qmcpack_4plyDisplay.Representation = 'Wireframe'
        qmcpack_4plyDisplay.AmbientColor = [1.0, 0.0, 0.0]
        qmcpack_4plyDisplay.ColorArrayName = [None, '']
        qmcpack_4plyDisplay.DiffuseColor = [1.0, 0.0, 0.0]
        qmcpack_4plyDisplay.LineWidth = 4.0
        qmcpack_4plyDisplay.SelectTCoordArray = 'None'
        qmcpack_4plyDisplay.SelectNormalArray = 'Normals'
        qmcpack_4plyDisplay.SelectTangentArray = 'None'
        qmcpack_4plyDisplay.OSPRayScaleArray = 'Normals'
        qmcpack_4plyDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
        qmcpack_4plyDisplay.SelectOrientationVectors = 'None'
        qmcpack_4plyDisplay.ScaleFactor = 0.30000000000000004
        qmcpack_4plyDisplay.SelectScaleArray = 'None'
        qmcpack_4plyDisplay.GlyphType = 'Arrow'
        qmcpack_4plyDisplay.GlyphTableIndexArray = 'None'
        qmcpack_4plyDisplay.GaussianRadius = 0.015
        qmcpack_4plyDisplay.SetScaleArray = ['POINTS', 'Normals']
        qmcpack_4plyDisplay.ScaleTransferFunction = 'PiecewiseFunction'
        qmcpack_4plyDisplay.OpacityArray = ['POINTS', 'Normals']
        qmcpack_4plyDisplay.OpacityTransferFunction = 'PiecewiseFunction'
        qmcpack_4plyDisplay.DataAxesGrid = 'GridAxesRepresentation'
        qmcpack_4plyDisplay.PolarAxes = 'PolarAxesRepresentation'
        qmcpack_4plyDisplay.SelectInputVectors = ['POINTS', 'Normals']
        qmcpack_4plyDisplay.WriteLog = ''

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        qmcpack_4plyDisplay.ScaleTransferFunction.Points = [-0.5773502588272095, 0.0, 0.5, 0.0, 0.5773502588272095, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        qmcpack_4plyDisplay.OpacityTransferFunction.Points = [-0.5773502588272095, 0.0, 0.5, 0.0, 0.5773502588272095, 1.0, 0.5, 0.0]


# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(tTKIcospheresFromPoints4)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')