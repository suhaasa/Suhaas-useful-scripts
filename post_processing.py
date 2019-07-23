# Ingrid Lan

# Automation script for post-processing my 3D simulations
# Makes use of both VMTK & Python VTK

# Reference for Python VTK scripting: https://conference.scipy.org/proceedings/scipy2015/pdfs/cory_quammen.pdf

import csv
import glob
import math
import os
import numpy as np
import time
import vtk
import vtk.util.numpy_support as nps

CONVERSION = 1333.2                 # conversion factor from mm Hg to dynes/cm^2

# Read a vtp file and return the polydata
def read_vtp(infile):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(infile)
    reader.Update()
    poly_data = reader.GetOutput()
    return poly_data


# Read a vtu file and return the vtu_data
def read_vtu(infile):
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(infile)
    vtu_reader.Update()
    vtu_data = vtu_reader.GetOutput()
    return vtu_data
  

# Sort cap VTP files into inflow / RPA branches / LPA branches. Obtain their names & cap areas.
def vtp_info():
    # Isolate just the .vtp's from all files in the mesh-surfaces directory
    filelist_raw = glob.glob(MESH_SURFACES_PATH)
    filelist_raw.sort()               # sort alphabetically
    filelist = []
    for trial in filelist_raw:
        if (trial[-4 : ] == ".vtp"):
          filelist.append(trial)

    # Sort caps into inflow / RPA branches / LPA branches. Store their names.
    rpa_names = []
    lpa_names = []
    inflow_name = ''

    inflow_tag = MESH_TAGS[0]
    rpa_branch_tag = MESH_TAGS[1]
    lpa_branch_tag = MESH_TAGS[2]

    for vtp_file in filelist:
        tail_name = vtp_file[len(MESH_SURFACES_PATH) - 1 : ]
        if (tail_name[ : len(rpa_branch_tag)] == rpa_branch_tag):
          rpa_names.append(vtp_file)
        elif (tail_name[ : len(lpa_branch_tag)] == lpa_branch_tag):
          lpa_names.append(vtp_file)
        elif (tail_name[ : len(inflow_tag)] == inflow_tag):
          inflow_name = vtp_file
    return {'rpa_names': rpa_names, 'lpa_names': lpa_names, 'inflow_name': inflow_name,
          'num_outlets': len(rpa_names) + len(lpa_names)}


# Return the normal slice at an interior point on the centerline
def interior_slice(interior_pt, norm_vec, all_results_vtu):
    # Use the centerline point & its normal vector to slice the .vtu model to return the perpendicular slice
    # at the centerline point.
    vtu_data = read_vtu(all_results_vtu)

    cutPlane = vtk.vtkPlane()
    cutPlane.SetOrigin(interior_pt)
    cutPlane.SetNormal(norm_vec)
    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(cutPlane)
    cutter.SetInputData(vtu_data)
    cutter.Update()

    # Use the connectivity filter to extract only regions connected to the interior point
    # This prevents slicing into other branches
    connect = vtk.vtkConnectivityFilter()
    connect.SetInputData(cutter.GetOutput())
    connect.SetExtractionModeToClosestPointRegion()
    connect.SetClosestPoint(interior_pt[0], interior_pt[1], interior_pt[2])
    connect.Update()

    vtu_slice = connect.GetOutput()
    return vtu_slice


# Write out a vtu slice - only used for testing!
# When loading the resulting vtu slice into Paraview, we're prompted to specify the reader to use -> choose VtkPolyData reader
def write_vtu_slice(vtu_slice, filename):
    slice_writer = vtk.vtkXMLPolyDataWriter()
    # slice_writer = vtk.vtkXMLUnstructuredGridWriter()
    slice_writer.SetInputData(vtu_slice)
    slice_writer.SetFileName(filename)
    slice_writer.Write()

# Integrate the vtu slice to obtain the flow & pressure over time at the given vtu_slice
def integrate_vtu_slice(vtu_slice, norm_vec, time_array, outfile):

    # Determine the area of each cell in the slice
    num_cells = vtu_slice.GetNumberOfCells()
    cell_areas = np.zeros(num_cells)
    for i in range(num_cells):
        cell = vtu_slice.GetCell(i)
        p0 = cell.GetPoints().GetPoint(0)
        p1 = cell.GetPoints().GetPoint(1)
        p2 = cell.GetPoints().GetPoint(2)
        cell_areas[i] = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
    total_area = sum(cell_areas)

    num_timepts = len(time_array)
    slice_flows = np.zeros((num_timepts, num_cells))
    slice_pressures = np.zeros(num_timepts)
    for i_time in range(num_timepts):
        velocity_name = 'velocity_' + time_array[i_time]
        pressure_name = 'pressure_' + time_array[i_time]
        
        velocity_array = vtu_slice.GetPointData().GetArray(velocity_name)
        pressure_array = vtu_slice.GetPointData().GetArray(pressure_name)

        for j_cell in range(num_cells):
            cell = vtu_slice.GetCell(j_cell)
            pts_cell = cell.GetPointIds()

            # For each cell, use trapezoidal integration to compute surface flow & pressure
            cell_flow = 0.0
            cell_pressure = 0.0

            for k_pt in xrange(0, pts_cell.GetNumberOfIds()):
                pt_id = pts_cell.GetId(k_pt)
                nodal_velocity = velocity_array.GetTuple3(pt_id)
                nodal_pressure = pressure_array.GetTuple1(pt_id)

                # Calculate the normal component of velocity (dot product)
                normal_velocity = sum(np.array(nodal_velocity) * norm_vec)
                                
                cell_flow += normal_velocity
                cell_pressure += nodal_pressure

            # Multiply the area by the mean normal velocity/mean pressure
            cell_flow *= cell_areas[j_cell] / float(pts_cell.GetNumberOfIds())
            cell_pressure *= cell_areas[j_cell] / float(pts_cell.GetNumberOfIds())

            # Store the cell's flow first, since the sign of all cell flows might need to be flipped if the
            # slice's normal vector points in the opposite direction as the flow.
            # Add the cell's pressure to the entire slice's pressure
            slice_flows[i_time, j_cell] = cell_flow
            slice_pressures[i_time] += cell_pressure
    
    # Divide slice pressures over time by the total area. Convert from dynes/cm2 to mm Hg
    slice_pressures = np.array(slice_pressures) / total_area / CONVERSION

    # The slice's normal vector could be pointing in the opposite direction as the flow.
    # If 5/10 of the slice flow's first time points are all negative, then take the negative of all flows
    if np.sum(np.sum(slice_flows, axis=1)[ : 10] < 0) > 5:
        slice_flows = -slice_flows

    # Sum up the contributon of all cells in the slice
    slice_flows = np.sum(slice_flows, axis=1)

    # Write out slice pressure (mm Hg) & flow at each time step of the last cardiac cycle
    string = "Time_Step Pressure_(mm_Hg) Flow_(mL/s)\n"
    for i in range(num_timepts):
        string += time_array[i] + " " + str(slice_pressures[i]) + " " + str(slice_flows[i]) + "\n"

    with open(outfile, 'w') as file:
        file.write(string)

    # Determine the period
    with open('../../pulmonary_targets.csv', 'r') as targets:
        targets_reader = csv.reader(targets)
        for row in targets_reader:
            if row[0] == 'Heart Rate':
                period = 60.0 / int(row[1])

    # Return the systolic/diastolic/mean pressures (mm Hg) & mean slice flow
    mean_slice_flow = np.sum(slice_flows) / period
    return (np.asarray([max(slice_pressures), min(slice_pressures), np.mean(slice_pressures)]), mean_slice_flow)


# Integrate a cap vtp to obtain the flow & pressure over time
def integrate_cap_vtp(results_vtp, cap_vtp, time_array, outfile):
    # First store all element IDs in the cap vtp
    cap_cell_ids = cap_vtp.GetCellData().GetArray('GlobalElementID')
    cap_cell_ids = nps.vtk_to_numpy(cap_cell_ids)

    # Determine the area of each cell in the cap vtp
    cap_num_cells = cap_vtp.GetNumberOfCells()
    cap_cell_areas = np.zeros(cap_num_cells)
    for i in range(cap_num_cells):
        cell = cap_vtp.GetCell(i)
        p0 = cell.GetPoints().GetPoint(0)
        p1 = cell.GetPoints().GetPoint(1)
        p2 = cell.GetPoints().GetPoint(2)
        cap_cell_areas[i] = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
    total_area = sum(cap_cell_areas)

    # Must loop through all cells in all_results.vtp (which contains solutions) and check
    # whether each cell is part of the cap vtp
    results_num_cells = results_vtp.GetNumberOfCells()
    results_cell_ids = results_vtp.GetCellData().GetArray('GlobalElementID')

    # Find the normal vector to the cap vtp
    normalGenerator = vtk.vtkPolyDataNormals()
    normalGenerator.SetInputData(cap_vtp)
    normalGenerator.ComputeCellNormalsOn()
    normalGenerator.ComputePointNormalsOff()
    normalGenerator.Update()
    normals_gen = normalGenerator.GetOutput().GetCellData().GetArray('Normals')
    normals = np.empty((0, 3))

    for i in range(cap_num_cells):
        norm_vec = np.asarray(normals_gen.GetTuple3(i))
        norm_vec_mag = np.sqrt(np.sum(norm_vec ** 2))
        normals = np.vstack((normals, norm_vec / norm_vec_mag))

    norm_vec = np.mean(normals, axis=0)

    num_timepts = len(time_array)
    cap_flows = np.zeros((num_timepts, cap_num_cells))
    cap_pressures = np.zeros(num_timepts)

    for i_time in range(num_timepts):
        velocity_name = 'velocity_' + time_array[i_time]
        pressure_name = 'pressure_' + time_array[i_time]
        
        velocity_array = results_vtp.GetPointData().GetArray(velocity_name)
        pressure_array = results_vtp.GetPointData().GetArray(pressure_name)

        for j_cell in range(results_num_cells):
            cell_id = results_cell_ids.GetTuple1(j_cell)

            # Only perform calculations if cell is part of the cap vtp
            idx = np.where(cell_id == cap_cell_ids)[0]
            if len(idx) > 0:
                cell = results_vtp.GetCell(j_cell)
                pts_cell = cell.GetPointIds()

                # For each cell, use trapezoidal integration to compute surface flow & pressure
                cell_flow = 0.0
                cell_pressure = 0.0

                for k_pt in xrange(0, pts_cell.GetNumberOfIds()):
                    pt_id = pts_cell.GetId(k_pt)
                    nodal_velocity = velocity_array.GetTuple3(pt_id)
                    nodal_pressure = pressure_array.GetTuple1(pt_id)

                    # Calculate the normal component of velocity (dot product)
                    normal_velocity = sum(np.array(nodal_velocity) * norm_vec)
                                    
                    cell_flow += normal_velocity
                    cell_pressure += nodal_pressure

                # Multiply the area by the mean normal velocity/mean pressure
                cell_flow *= cap_cell_areas[idx[0]] / float(pts_cell.GetNumberOfIds())
                cell_pressure *= cap_cell_areas[idx[0]] / float(pts_cell.GetNumberOfIds())

                # Store the cell's flow first, since the sign of all cell flows might need to be flipped if the
                # slice's normal vector points in the opposite direction as the flow.
                # Add the cell's pressure to the entire slice's pressure
                cap_flows[i_time, idx[0]] = cell_flow
                cap_pressures[i_time] += cell_pressure
        
    # Divide cap pressures over time by the total area. Convert from dynes/cm2 to mm Hg
    cap_pressures = np.array(cap_pressures) / total_area / CONVERSION

    # The slice's normal vector could be pointing in the opposite direction as the flow.
    # If 5/10 of the slice flow's first time points are all negative, then take the negative of all flows
    if np.sum(np.sum(cap_flows, axis=1)[ : 10] < 0) > 5:
        cap_flows = -cap_flows

    # Sum up the contributon of all cells in the slice
    cap_flows = np.sum(cap_flows, axis=1)

    # Write out cap pressure (mm Hg) & flow at each time step of the last cardiac cycle
    string = "Time_Step Pressure_(mm_Hg) Flow_(mL/s)\n"
    for i in range(num_timepts):
        string += time_array[i] + " " + str(cap_pressures[i]) + " " + str(cap_flows[i]) + "\n"

    if outfile != '':
        with open(outfile, 'w') as file:
            file.write(string)

    # Return the systolic/diastolic/mean pressures (mm Hg)
    return np.asarray([max(cap_pressures), min(cap_pressures), np.mean(cap_pressures)])


# Post-process the all_results.vtu and all_results.vtp files
def post_process(results_vtu, results_vtp):
    inflow_tag = MESH_TAGS[0]
    inflow_vtp = MESH_SURFACES_PATH[ : -1] + inflow_tag + '.vtp'
    results_vtp = read_vtp(results_vtp)
    inflow_vtp = read_vtp(inflow_vtp)

    # Read in the RPA, LPA slice coordinates & norm vectors from 'slice_locs.txt'
    slice_coords = np.zeros((2, 3))             # RPA/LPA x XYZ
    norm_vecs = np.zeros((2, 3))
    with open('../../slice_locs.txt', 'r') as f:
        lines = f.read().splitlines()
        for i in range(len(lines)):
            line_split = lines[i].split()
            row = int(math.floor(i / 2))
            if i % 2 == 0:                      # even lines correspond to slice coordinates
                slice_coords[row, :] = np.asarray([float(line_split[1]), float(line_split[2]), float(line_split[3])])
            else:                               # odd lines correspond to unit normal vectors
                norm_vecs[row, :] = np.asarray([float(line_split[1]), float(line_split[2]), float(line_split[3])])

    # Determine the number of time points with solutions
    time_array = []
    for i in range(0, results_vtp.GetPointData().GetNumberOfArrays()):
        array_name = results_vtp.GetPointData().GetArray(i).GetName()
        if array_name[0 : 8] == 'velocity':
            time_array.append(array_name[9 : ])
      
    num_timepts = len(time_array)
    print('Number of timepoints: ' + str(num_timepts))

    # Generate RPA, LPA slices from all_results_vtu. Then integrate each slice to determine the
    # slice pressure & flow at each time step in the last cardiac cycle
    filenames = ['rpa_slice_results.txt', 'lpa_slice_results.txt']
    sim_pressures = np.zeros((len(filenames), 3))            # RPA/LPA x systolic/diastolic/mean pressure  
    mean_slice_flows = np.zeros(len(filenames))
    mean_slice_pressures = np.zeros(len(filenames))     
    for i in range(len(filenames)):
        vtu_slice = interior_slice(slice_coords[i, :], norm_vecs[i, :], results_vtu)
        slice_pressures, slice_flow = integrate_vtu_slice(vtu_slice, norm_vecs[i, :], time_array, filenames[i])
        sim_pressures[i, :] = slice_pressures
        mean_slice_flows[i] = slice_flow
        mean_slice_pressures[i] = slice_pressures[-1]

    # Integrate inflow vtp to determine the MPA pressure & flow at each time step in the last cardiac cycle
    mpa_pressures = integrate_cap_vtp(results_vtp, inflow_vtp, time_array, 'inflow_results.txt')
    sim_pressures = np.vstack((mpa_pressures, sim_pressures))

    # Write out to 'sim_pressures.txt'
    string = ''
    for i in range(sim_pressures.shape[0]):
        for j in range(sim_pressures.shape[1]):
            string += '{:.4f}'.format(sim_pressures[i, j]) + ' '
        string += '\n'
    with open('sim_pressures.txt', 'w+') as sim_results:
        sim_results.write(string)

    # Compute the mean pressure at each RPA & LPA outlet
    caps = vtp_info()
    rpa_outlets = caps['rpa_names']
    rpa_outlet_pressures = np.zeros(len(caps['rpa_names']))
    for i in range(len(rpa_outlets)):
        cap_vtp = read_vtp(rpa_outlets[i])
        
        # Integrate inflow vtp to determine the pressure & flow at each time step in the last cardiac cycle
        rpa_outlet_pressures[i] = integrate_cap_vtp(results_vtp, cap_vtp, time_array, '')[-1]

    lpa_outlets = caps['lpa_names']
    lpa_outlet_pressures = np.zeros(len(caps['lpa_names']))
    for i in range(len(lpa_outlets)):
        cap_vtp = read_vtp(lpa_outlets[i])
        
        # Integrate inflow vtp to determine the pressure & flow at each time step in the last cardiac cycle
        lpa_outlet_pressures[i] = integrate_cap_vtp(results_vtp, cap_vtp, time_array, '')[-1]

    # Determine resistances to represent the distal vasculature past the proximal stenoses
    R_3d_1 = CONVERSION * (mean_slice_pressures[0] - np.mean(rpa_outlet_pressures)) / mean_slice_flows[0]
    R_3d_2 = CONVERSION * (mean_slice_pressures[1] - np.mean(lpa_outlet_pressures)) / mean_slice_flows[1]
    if R_3d_1 < 0:
        R_3d_1 = 0
    if R_3d_2 < 0:
        R_3d_2 = 0
    return (R_3d_1, R_3d_2)


if __name__ == "__main__":
    print("** Post-processing deformable simulation results...\n")
    with open('postprocessing_command.txt', 'r') as command:
        lines = command.read().splitlines()
        command_string = lines[0]
        tuning_iter = int(lines[1])
    print(command_string)
    os.system(command_string)

    # Will be in n-procs_case folder
    while not os.path.isfile('all_results.vtp'):
        time.sleep(60)   

    command_string = 'cd ..'
    print(command_string)
    os.chdir('..')

    # Move deformable results into a new directory
    results_dir = 'deformable_results_' + str(tuning_iter)
    print('\nMoving deformable results into ' + results_dir)
    command_string = 'mkdir ' + results_dir
    print(command_string)
    os.system(command_string)

    procs_folder = glob.glob('*-procs_*')[0]
    output_file = glob.glob('*.o*')[0]
    error_file = glob.glob('*.e*')[0]
    flow_file = glob.glob('*.flow')[0]

    command_string = 'mv varwallprop.vtp displacement.vtp rcrt.dat bct.vtp bct.dat restart.0.1 geombc.dat.1 ' + procs_folder + ' ' + \
                     output_file + ' ' + error_file + ' ' + flow_file + ' ' + results_dir
    print(command_string)
    os.system(command_string)

    command_string = 'cd ' + results_dir
    print(command_string)
    os.chdir(results_dir)

    # Determine the systolic/mean/diastolic pressures of post-stenosis RPA & LPA slices
    procs_folder = glob.glob('*-procs_*')[0]
    results_vtu = procs_folder + '/all_results.vtu'
    results_vtp = procs_folder + '/all_results.vtp'

    # Read in path to all mesh surfaces
    MESH_SURFACES_PATH = ''
    with open('../../mesh_surf_path.txt', 'r') as infile:
        MESH_SURFACES_PATH = infile.readline()

    MESH_TAGS = []
    with open('../../mesh_tags.txt', 'r') as infile:
        lines = infile.read().splitlines()
        MESH_TAGS = lines[0].split()
    
    # Post-process all_results.vtu and all_results.vtp
    R_3d_1, R_3d_2 = post_process(results_vtu, results_vtp)

    with open('distal_resistances.txt', 'w+') as file:
        file.write(str(R_3d_1) + ' ' + str(R_3d_2))

    print('Post-processing completed')
