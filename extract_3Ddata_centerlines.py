##wgyang 2019.3
##This python script extracts 3D data in the 3D simulation based on groups divided by vmtk. It is used for comparing 1D results
##
import math
import sys
import os
from os import path, makedirs
#import pathlib
import numpy as np
import numpy.linalg as la

import vtk
from vmtk import vtkvmtk,vmtkscripts
from vtk.util import numpy_support
import vtk.util.numpy_support as nps
import os
import numpy as np
#pv_path="/home/jeff/ParaView-5.6.0-MPI-Linux-64bit/lib"
#sys.path.append('%s/python2.7/site-packages' % pv_path)
#sys.path.append('%s/python2.7/site-packages/vtkmodules' % pv_path)
#sys.path.append(pv_path)
#from paraview import simple as pvs
############################################################################################################
## user inputs
ModelName="2017_398_def_PCMRI"
centerline_file="PH_2017_398_cl.vtp"
averageWSS_file="average_result.vtp"
tstart=6000
tfin=8000
incr=40
tcyle=0.76
#0: segment inlet, 0.5: middle, 1:segment outlet  
seg_location=0.5
iflowpress=1
igenarea=0
igenWSS=0
#output file names
output_flowfile=ModelName+"_flow.dat"
output_pressfile=ModelName+"_press.dat"
output_areafile=ModelName+"_area.dat"
output_WSSfile=ModelName+"_WSS.dat"

###########################################################################################################


# Read a vtp file and return the polydata
def read_polydata(filename, datatype=None):
    """
    Load the given file, and return a vtkPolyData object for it.
    Args:
        filename (str): Path to input file.
        datatype (str): Additional parameter for vtkIdList objects.
    Returns:
        polyData (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Output data.
    """

    # Check if file exists
    if not path.exists(filename):
        raise RuntimeError("Could not find file: %s" % filename)

    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get reader
    if fileType == 'stl':
        reader = vtk.vtkSTLReader()
        reader.MergingOn()
    elif fileType == 'vtk':
        reader = vtk.vtkPolyDataReader()
    elif fileType == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif fileType == 'vts':
        reader = vtk.vtkXMinkorporereLStructuredGridReader()
    elif fileType == 'vtr':
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fileType == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fileType == "vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fileType == "np" and datatype == "vtkIdList":
        result = np.load(filename).astype(np.int)
        id_list = vtk.vtkIdList()
        id_list.SetNumberOfIds(result.shape[0])
        for i in range(result.shape[0]):
            id_list.SetId(i, result[i])
        return id_list
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Read
    reader.SetFileName(filename)
    reader.Update()
    polydata = reader.GetOutput()

    return polydata

def write_polydata(input_data, filename, datatype=None):
    """
    Write the given input data based on the file name extension.
    Args:
        input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Input data.
        filename (str): Save path location.
        datatype (str): Additional parameter for vtkIdList objects.
    """
    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get writer
    if fileType == 'stl':
        writer = vtk.vtkSTLWriter()
    elif fileType == 'vtk':
        writer = vtk.vtkPolyDataWriter()
    elif fileType == 'vts':
        writer = vtk.vtkXMLStructuredGridWriter()
    elif fileType == 'vtr':
        writer = vtk.vtkXMLRectilinearGridWriter()
    elif fileType == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif fileType == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif fileType == "vti":
        writer = vtk.vtkXMLImageDataWriter()
    elif fileType == "np" and datatype == "vtkIdList":
        output_data = np.zeros(input_data.GetNumberOfIds())
        for i in range(input_data.GetNumberOfIds()):
            output_data[i] = input_data.GetId(i)
        output_data.dump(filename)
        return
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Set filename and input
    writer.SetFileName(filename)
    writer.SetInputData(input_data)
    writer.Update()

    # Write
    writer.Write()


# Calculate the centroid of a vtp file
def centroid(infile):
    poly_data = read_polydata(infile)
    x_list = []
    y_list = []
    z_list = []
    for i in range(poly_data.GetNumberOfPoints()):
        x_list.append(float(poly_data.GetPoints().GetPoint(i)[0]))
        y_list.append(float(poly_data.GetPoints().GetPoint(i)[1]))
        z_list.append(float(poly_data.GetPoints().GetPoint(i)[2]))

    return [np.mean(x_list), np.mean(y_list), np.mean(z_list)]

################################################################################
centerline_data=read_polydata(centerline_file)
#calculate Frenet tangential vector
centergeometry=vmtkscripts.vmtkCenterlineGeometry()
centergeometry.Centerlines=centerline_data
centergeometry.Execute()
centerline_data=centergeometry.Centerlines
write_polydata(centerline_data,centerline_file)



num_pts=centerline_data.GetNumberOfPoints()
print "Number of Points:", centerline_data.GetNumberOfPoints()

print "Number of arrays:", centerline_data.GetCellData().GetNumberOfArrays()

num_cells=centerline_data.GetNumberOfCells()
print "Number of Cells:", centerline_data.GetNumberOfCells()

###read cell data, for each cell (line), the lists record its centerline id (n lines starting from the inlet to outlets), blanking (0 non bifurcation, 1 bifurcation), 
###group ids (vessels are splitted into single segments/branches and a bifurcation region, 
###if a vessel shared by multiple centerlines, there multiple elements sharing the same group id
###if a segment is the terminal segment without any child elements, its group id is unique.
###for a bifurcation region, the elemens' blanking is 1.

celldata=centerline_data.GetCellData().GetArray("CenterlineIds")
centerline_list = nps.vtk_to_numpy(celldata)

celldata=centerline_data.GetCellData().GetArray("Blanking")
blank_list = nps.vtk_to_numpy(celldata)

celldata=centerline_data.GetCellData().GetArray("GroupIds")
group_list = nps.vtk_to_numpy(celldata)

celldata=centerline_data.GetCellData().GetArray("TractIds")
tract_list = nps.vtk_to_numpy(celldata)


num_path=centerline_list[-1]+1
print "number of paths=", num_path

num_group=max(group_list)+1
print "number of groups=",num_group


###path_elems[i] records the element(line) indices for centerline id=i
path_elems=[]

for i in range(0,num_path):
 path_elems.append([])

for i in range(0,num_path):
  for j in range(0,num_cells):
    if i==centerline_list[j]: 
     path_elems[i].append(j)

for i in range(0,num_path):
 print "centerline",i,"groups ids",path_elems[i]

###group_elems[i] records the element(line) indices for group id=i
group_elems=[]
for i in range(0,num_group):
  group_elems.append([])


for i in range(0,num_cells):
#  print "cellid=",i
#  print "groupid=",group_list[i]
  group_elems[group_list[i]].append(i)

for i in range(0,num_group):
 print "group",i, "element ids", group_elems[i]

###path_elems[i] records the element(line) indices for centerline id=i
path_elems=[]

for i in range(0,num_path):
 path_elems.append([])

for i in range(0,num_path):
  for j in range(0,num_cells):
    if i==centerline_list[j]: 
     path_elems[i].append(j)

for i in range(0,num_path):
 print "centerline",i,"groups ids",path_elems[i]

###group_elems[i] records the element(line) indices for group id=i
group_elems=[]
for i in range(0,num_group):
  group_elems.append([])


for i in range(0,num_cells):
#  print "cellid=",i
#  print "groupid=",group_list[i]
  group_elems[group_list[i]].append(i)


group_terminal=[0]*num_group
num_outlet=0
num_bif=0
tmp=len(group_elems[0])
for i in range(0,num_group):
 #this is a bifurcation
 if blank_list[group_elems[i][0]]==1:
  group_terminal[i]=2
  num_bif=num_bif+1
 # this is a terminal segment
 if len(group_elems[i])==1:
  group_terminal[i]=1
  num_outlet=num_outlet+1
 if len(group_elems[i])>tmp and blank_list[group_elems[i][0]]!=1:
  tmp=len(group_elems[i])

print "group_terminal=",group_terminal


if tmp!=len(group_elems[0]) or tmp!=num_outlet or num_path!=num_outlet:
 print "warning: inlet group id is not 0 or number of centerlines is not equal to the number of outlets"
 exit()


pointdata=centerline_data.GetPointData().GetArray("MaximumInscribedSphereRadius")
points_maxR=nps.vtk_to_numpy(pointdata)

pointdata=centerline_data.GetPointData().GetArray("FrenetTangent")
points_tangent=nps.vtk_to_numpy(pointdata)

#For each group, the center of roi is recorded
group_roi_center=[]
group_roi_maxR=[]
group_roi_tan=[]
for i in range(0, num_group):
   ids=vtk.vtkIdList()
   centerline_data.GetCellPoints(group_elems[i][0],ids)
   num_ids=ids.GetNumberOfIds()
   print "group=",i,"num of points=",num_ids
   path_dist=np.zeros(num_ids)
   for pathptid in range(1, num_ids):
      id1=ids.GetId(pathptid)
      id2=ids.GetId(pathptid-1)
      pt1=np.array(centerline_data.GetPoints().GetPoint(id1))
      pt2=np.array(centerline_data.GetPoints().GetPoint(id2))
      path_dist[pathptid]=path_dist[pathptid-1]+np.linalg.norm(pt1-pt2)
   
   dist2target=abs(path_dist-path_dist[num_ids-1]*seg_location) 
   index=np.argmin(dist2target)
  # index=int(round((num_ids-1)*seg_location))
   tmpid=ids.GetId(index)
   
   
   tmpR=points_maxR[tmpid]
   pt1=np.array(centerline_data.GetPoints().GetPoint(tmpid))
   
   group_roi_center.append(pt1)
   group_roi_maxR.append(tmpR)
   if np.linalg.norm(points_tangent[tmpid])<1e-6:
     
     if index<num_ids-1:
       pt2=np.array(centerline_data.GetPoints().GetPoint(ids.GetId(index+1)))
           
     else:
       pt2=np.array(centerline_data.GetPoints().GetPoint(ids.GetId(index-1)))
     dx=np.linalg.norm(pt2-pt1)
     tmptan=(pt2-pt1)/dx
     print "tangent finite diff",tmptan
     group_roi_tan.append(tmptan)
   else:
     group_roi_tan.append(points_tangent[tmpid])

#print "group_roi_center:",group_roi_center
#print "group_roi_maxR:",group_roi_maxR
#print "group_roi_tan:",group_roi_tan


if iflowpress==1:
 group_flow=[]
 group_pressure=[]


 for i in range(0,num_group):
   group_flow.append([])
   group_pressure.append([])



 for i in range(tstart,tfin+incr,incr):
   filename=ModelName+ "_%05d.vtu"%i 
   vtu_data=read_polydata(filename)
  
   for j in range(0,num_group):
    if group_terminal[j]!=2:
         
      plane = vtk.vtkPlane()
      plane.SetOrigin(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
      plane.SetNormal(group_roi_tan[j][0], group_roi_tan[j][1], group_roi_tan[j][2])

      cutter = vtk.vtkCutter()
      cutter.SetCutFunction(plane)
      cutter.SetInputData(vtu_data)
      cutter.Update()
      vtu_slice=cutter.GetOutput()
     
    #  write_polydata(vtu_slice,"group_slice"+str(j)+".vtp")

      sphere = vtk.vtkSphereSource()
      sphere.SetCenter(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
      sphere.SetRadius(group_roi_maxR[j]*1.5)
      sphere.Update()
     
      implicitPolyDataDistance =vtk.vtkImplicitPolyDataDistance()
      implicitPolyDataDistance.SetInput(sphere.GetOutput())
     
      signedDistances=vtk.vtkFloatArray()
      signedDistances.SetNumberOfComponents(1)
      signedDistances.SetName("SignedDistances")
     # print"cutter num points",vtu_slice.GetNumberOfPoints()
      if vtu_slice.GetNumberOfPoints()==0:
          print "0 sliced point for group=",j
          exit()
      for pointId in range(vtu_slice.GetNumberOfPoints()):
          p=vtu_slice.GetPoint(pointId)
          signedDistance=implicitPolyDataDistance.EvaluateFunction(p)
          signedDistances.InsertNextValue(signedDistance)
     
      # default scalar field pressure is replaced by setscalars(signeddistances), get pressure array first  
      pressure_array = vtu_slice.GetPointData().GetArray("pressure")
      vtu_slice.GetPointData().SetScalars(signedDistances) 
      vtu_slice.GetPointData().AddArray(pressure_array)
      clipper=vtk.vtkClipDataSet()
      clipper.SetInputData(vtu_slice)
      clipper.InsideOutOn()
      clipper.SetValue(0.0)
      clipper.Update()
    

      vtu_slice=clipper.GetOutput()
      print "number of points clip",vtu_slice.GetNumberOfPoints()
      if vtu_slice.GetNumberOfPoints()==0:
         print "0 clipped point for group=",j
         exit()

      vtu_slice_file="group"+str(j)+"_slice"+".vtu"
      print "t=",i,vtu_slice_file
     
      vn=[0.0,0.0,0.0]
      slice_area=0.0
      slice_flow=0.0
      slice_pressure=0.0
      cell = vtu_slice.GetCell(0)
      p0 = cell.GetPoints().GetPoint(0)
      p1 = cell.GetPoints().GetPoint(1)
      p2 = cell.GetPoints().GetPoint(2)
      vtk.vtkTriangle().ComputeNormal(p0,p1,p2,vn)
      pressure_array = vtu_slice.GetPointData().GetArray("pressure")
      velocity_array = vtu_slice.GetPointData().GetArray("velocity")
      for cellId in range(vtu_slice.GetNumberOfCells()):
         cell = vtu_slice.GetCell(cellId)
         p0 = cell.GetPoints().GetPoint(0)
         p1 = cell.GetPoints().GetPoint(1)
         p2 = cell.GetPoints().GetPoint(2) 
         cell_area=vtk.vtkTriangle().TriangleArea(p0, p1, p2)
         slice_area=slice_area+cell_area
         nodeids = cell.GetPointIds()
         v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
         v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
         v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))  
         slice_flow=slice_flow+sum((v1+v2+v3)/3.0*vn)*cell_area # 1 point Gauss quad rule
         p0 = np.array(pressure_array.GetTuple(nodeids.GetId(0)))
         p1 = np.array(pressure_array.GetTuple(nodeids.GetId(1)))
         p2 = np.array(pressure_array.GetTuple(nodeids.GetId(2)))
         slice_pressure=slice_pressure+np.mean([p0,p1,p2])*cell_area
        
     
        
      slice_pressure=slice_pressure/slice_area
      print "slice_area=",slice_area,"flow",slice_flow,"pressure",slice_pressure
      group_flow[j].extend([slice_flow])
      group_pressure[j].extend([slice_pressure])
      write_polydata(vtu_slice,vtu_slice_file)
    else:
      # skip bifurcation group
      group_flow[j].extend([0.0])
      group_pressure[j].extend([0.0])
      

 flowfile = open(output_flowfile, "w")
 pressfile =open(output_pressfile,"w")
 for i in range(0,num_group):
   for j in range(0,len(group_flow[0])):
    flowfile.write(str(group_flow[i][j])+" ")
    pressfile.write(str(group_pressure[i][j])+" ")
   flowfile.write("\n")
   pressfile.write("\n")
 flowfile.close()
 pressfile.close()


###process area WSS 
if igenarea==1:
 group_area=[]
 for i in range(0,num_group):
   group_area.append([])
 print "process displacement data"
 for i in range(tstart,tfin+incr,incr):
   filename=ModelName+ "_%05d.vtp"%i 
   filename2=ModelName+ "_%05d_realdisp.vtp"%i 
   vtp_data=read_polydata(filename) 
   if i==tstart:
     reference_disp=vtk.vtkDoubleArray()
     reference_disp=vtp_data.GetPointData().GetArray("displacement")
   
   current_disp=vtk.vtkDoubleArray()
   current_disp=vtp_data.GetPointData().GetArray("displacement") 

 
   numPts=vtp_data.GetNumberOfPoints()
   realdisplacement=vtk.vtkDoubleArray()
   realdisplacement.SetNumberOfComponents(3)
   realdisplacement.Allocate(numPts,10000)
   realdisplacement.SetNumberOfTuples(numPts)
   realdisplacement.SetName("realdisplacement")

   for ptid in xrange(0,numPts):
     ref_pt_disp=reference_disp.GetTuple3(ptid)
     cur_pt_disp=current_disp.GetTuple3(ptid)
     realdisplacement.SetTuple3(ptid,cur_pt_disp[0]-ref_pt_disp[0],cur_pt_disp[1]-ref_pt_disp[1],cur_pt_disp[2]-ref_pt_disp[2]) 
     vtp_data.GetPointData().AddArray(realdisplacement) 
     pt=np.array(vtp_data.GetPoints().GetPoint(ptid))
     newpt=pt+(np.array(cur_pt_disp)-np.array(ref_pt_disp))
     vtp_data.GetPoints().SetPoint(ptid,newpt)
   write_polydata(vtp_data,filename2)



 for i in range(tstart,tfin+incr,incr):
   filename2=ModelName+ "_%05d_realdisp.vtp"%i 
   vtp_data=read_polydata(filename2)  
    
   for j in range(0,num_group):
    if group_terminal[j]!=2:
      plane = vtk.vtkPlane()
      plane.SetOrigin(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
      plane.SetNormal(group_roi_tan[j][0], group_roi_tan[j][1], group_roi_tan[j][2])

      cutter = vtk.vtkCutter()
      cutter.SetCutFunction(plane)
      cutter.SetInputData(vtp_data)
      cutter.Update()
      vtp_slice=cutter.GetOutput()
    #  write_polydata(vtp_slice,"group"+str(j)+"_slice_realdisp.vtp")
      sphere = vtk.vtkSphereSource()
      sphere.SetCenter(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
      sphere.SetRadius(group_roi_maxR[j]*2.0)
      sphere.Update()
    #  print"group=",j,"sphere center=",group_roi_center[j][0],group_roi_center[j][1],group_roi_center[j][2],"radius=",group_roi_maxR[j]*1.5,"group_roi_maxR=",group_roi_maxR[j]
      implicitPolyDataDistance =vtk.vtkImplicitPolyDataDistance()
      implicitPolyDataDistance.SetInput(sphere.GetOutput())
     
      signedDistances=vtk.vtkFloatArray()
      signedDistances.SetNumberOfComponents(1)
      signedDistances.SetName("SignedDistances")
      print"displacement area cutter num points",vtp_slice.GetNumberOfPoints()
      if vtp_slice.GetNumberOfPoints()==0:
          exit()
      for pointId in range(vtp_slice.GetNumberOfPoints()):
          p=vtp_slice.GetPoint(pointId)
          signedDistance=implicitPolyDataDistance.EvaluateFunction(p)
          signedDistances.InsertNextValue(signedDistance)
     
      vtp_slice.GetPointData().SetScalars(signedDistances) 
      clipper=vtk.vtkClipDataSet()
      clipper.SetInputData(vtp_slice)
      clipper.InsideOutOn()
      clipper.SetValue(0.0)
      clipper.Update()
      vtu_clip=clipper.GetOutput()
      print "displacement number of points clip",vtu_clip.GetNumberOfPoints()

      if vtu_clip.GetNumberOfPoints()==0:
         exit()

      vtu_clip_file="group"+str(j)+"_slice_realdisp"+".vtu" 
      write_polydata(vtu_clip,vtu_clip_file)
     
      centerpt=np.array([0.0,0.0,0.0])
      for ptid in range(0,vtu_clip.GetNumberOfPoints()):
         centerpt=centerpt+np.array(vtu_clip.GetPoints().GetPoint(ptid))
        
      centerpt=centerpt/vtu_clip.GetNumberOfPoints()
     
      slice_area=0.0
     
     # print "clip ring num cell=",vtu_clip.GetNumberOfCells()
      for cellid in range(0,vtu_clip.GetNumberOfCells()):
         cell = vtu_clip.GetCell(cellid)
         p0 = np.array(cell.GetPoints().GetPoint(0))
         p1 = np.array(cell.GetPoints().GetPoint(1))
   
         slice_area =slice_area+ vtk.vtkTriangle().TriangleArea(p0, p1, centerpt)

      print "t=",i," group",j, "cut area=",slice_area
      group_area[j].extend([slice_area])
    else:
      group_area[j].extend([0.0])


 areafile = open(output_areafile, "w")
 for i in range(0,num_group):
   for j in range(0,len(group_area[0])):
    areafile.write(str(group_area[i][j])+" ")
   areafile.write("\n")
 areafile.close()




if igenWSS==1:
 group_WSS=np.zeros(num_group)
 group_area=np.zeros(num_group)
 print "process WSS data"
 print"split branch and calculate WSS for each group"
 modeldata=read_polydata(averageWSS_file)


  
 for i in range(0,num_group):
  if group_terminal[i]!=2: 
    clipgroup = vmtkscripts.vmtkBranchClipper()
    clipgroup.Surface=modeldata
    clipgroup.Centerlines=centerline_data
    clipgroup.GroupIds=[i]
    clipgroup.Execute()
    clipgroup_output = clipgroup.Surface
    write_polydata(clipgroup_output,"group_"+str(i)+"_avg.vtp")

    num_pts=clipgroup_output.GetPoints().GetNumberOfPoints()
    num_cells=clipgroup_output.GetNumberOfCells()
    print"clip group=",i,"num_pts=",num_pts,"num_cells=",num_cells
    groupid_array=nps.vtk_to_numpy(clipgroup_output.GetPointData().GetArray("GroupIds"))

    WSS_array=nps.vtk_to_numpy(clipgroup_output.GetPointData().GetArray("vTAWSS_wss")) 
 
 
    for cellId in range(0,num_cells):

      cell =clipgroup_output.GetCell(cellId)
      p0 = cell.GetPoints().GetPoint(0)
      p1 = cell.GetPoints().GetPoint(1)
      p2 = cell.GetPoints().GetPoint(2) 
      nodeids = cell.GetPointIds()     
      cell_area=vtk.vtkTriangle().TriangleArea(p0, p1, p2)
      
      v1 = np.array(WSS_array[nodeids.GetId(0)])
      v2 = np.array(WSS_array[nodeids.GetId(1)])
      v3 = np.array(WSS_array[nodeids.GetId(2)]) 
      if abs(np.mean([v1,v2,v3]))>1e-7: 
         group_WSS[i]=group_WSS[i]+np.mean([v1,v2,v3])*cell_area # 1 point Gauss quad rule
         group_area[i]=group_area[i]+cell_area

    print"group area=",group_area[i]  
 
 for i in range(0,num_group):
   if group_area[i]>0.0:
       group_WSS[i]=group_WSS[i]/group_area[i]     

 
 WSSfile = open(output_WSSfile, "w")
 for i in range(0,num_group):
    WSSfile.write(str(group_WSS[i])+" "+str(group_area[i])+"\n")
 WSSfile.close()

