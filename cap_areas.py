import glob
import vtk
from os import path, makedirs

#--------------------Define mesh information-------------------------#
MESH_SURFACES_PATH = '/home/suhaasa/Suhaas/adult_heart_2_mgs/Meshes/full_model_lca95_3col_28um-mesh-complete/mesh-surfaces/*'
FILENAME = 'cap_areas_adult_heart_2.csv'
inflow_tag = 'ao'
rpa_branch_tag = 'RCA'
lpa_branch_tag = 'L'
rsa_branch_tag = 'RSA'


def find_vtp_area(infile):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(infile)
    reader.Update()
    poly = reader.GetOutputPort()
    masser = vtk.vtkMassProperties()
    masser.SetInputConnection(poly)
    masser.Update()
    return masser.GetSurfaceArea()


def vtp_info():
    # Isolate just the .vtp's from all files in the mesh-surfaces directory
    filelist_raw = glob.glob(MESH_SURFACES_PATH)
    filelist_raw.sort()               # sort alphabetically
    filelist = []
    for trial in filelist_raw:
        if (trial[-4 : ] == ".vtp"):
          filelist.append(trial)

    # Sort caps into inflow / RPA branches / LPA branches. Store their names & cap areas.
    rpa_areas = []
    rpa_names = []
    lpa_areas = []
    lpa_names = []
    rsa_names = []
    rsa_areas = []
    inflow_area = []
    inflow_name = []

    for vtp_file in filelist:
        tail_name = vtp_file[len(MESH_SURFACES_PATH) - 1 : ]
        if (tail_name[ : len(rpa_branch_tag)] == rpa_branch_tag):
          rpa_areas.append(find_vtp_area(vtp_file))
          rpa_names.append(vtp_file)
        elif (tail_name[ : len(lpa_branch_tag)] == lpa_branch_tag):
          lpa_areas.append(find_vtp_area(vtp_file))
          lpa_names.append(vtp_file)
        elif (tail_name[ : len(rsa_branch_tag)] == rsa_branch_tag):
          rsa_areas.append(find_vtp_area(vtp_file))
          rsa_names.append(vtp_file)
        elif (tail_name[ : len(inflow_tag)] == inflow_tag):
          inflow_area.append(find_vtp_area(vtp_file))
          inflow_name.append(vtp_file)
    return {'rca_areas': rpa_areas, 'rca_names': rpa_names, 'lca_areas': lpa_areas, 
          'lca_names': lpa_names, 'rsa_areas': rsa_areas, 'rsa_names': rsa_names, 'inflow_area': inflow_area, 'inflow_name': inflow_name,
          'num_outlets': len(rpa_areas) + len(lpa_areas)}

def print_to_file(caps,filename):
    rca_areas = caps.get('rca_areas')
    rca_names = caps.get('rca_names')
    lca_areas = caps.get('lca_areas')
    lca_names = caps.get('lca_names')
    rsa_areas = caps.get('rsa_areas')
    rsa_names = caps.get('rsa_names')
    inflow_area = caps.get('inflow_area')
    inflow_name = caps.get('inflow_name')
    
    outfile = open(filename,'w')
    out_string = 'Cap,Area' + '\n'
    outfile.write(out_string)
    for i in range(0,len(lca_names)):
      out_string = lca_names[i] + ',' + str(lca_areas[i]) + '\n'
      outfile.write(out_string)
    for i in range(0,len(rca_names)):
      out_string = rca_names[i] + ',' + str(rca_areas[i]) + '\n'
      outfile.write(out_string)
    for i in range(0,len(rsa_names)):
      out_string = rsa_names[i] + ',' + str(rsa_areas[i]) + '\n'
      outfile.write(out_string)
    for i in range(0,len(inflow_name)):
      out_string = inflow_name[i] + ',' + str(inflow_area[i]) + '\n'
      outfile.write(out_string)
    outfile.close()
    print('wrote file.')


if __name__ == "__main__":

    # Sort cap VTP files into inflow / RPA branches / LPA branches. Obtain their names & cap areas.
    caps = vtp_info()
    print_to_file(caps,FILENAME)
