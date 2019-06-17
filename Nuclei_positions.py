#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

#--------------------------------------------------------------------------------------------------
# Should this be implemented in Fortran?Yes Fortran for PF and C for CA 
# I would need the post-processing steps like Gradient and rotation calculations
# But it would be better if done using the HDF5 data
# For that it is better to do in Python, it seems
# All of this should be done with restart method in DAMASK (need to ask ChuanLai how to do in DAMASK)
# the method based on restart is only possible for single processor based outputs
# Or work on getting CA data to DAMASK (then the same calculation can be done in CA)
#--------------------------------------------------------------------------------------------------

import numpy as np
import h5py
import damask
from dadf5 import DADF5
import dadf5
import argparse
from scipy import spatial

#--------------------------------------------------------------------------------------------------
# File to be read and modified with calculations 
#--------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Creating histograms')
parser.add_argument('file',nargs='+',help='filename where the data lies')
parser.add_argument('--geom',nargs='+',help='geom file name, required because HDF5 file doesnt have this info')
parser.add_argument('--restart',nargs='+',help='restart file for next DAMASK simulation')
args = parser.parse_args() #filenames
print(args.file)
print(args.geom)
print(args.restart)

if args.restart:
    restart_file = h5py.File(args.restart[0],'a')
#--------------------------------------------------------------------------------------------------
# Use dadf5 to get dislocation density gradients 
# Grain_rotation would already be added in the HDF5 file
# Add calculation for calculating 'possible' subgrain size at each material point
#--------------------------------------------------------------------------------------------------
#constants for calculations: (should be changed by the user)
G = 33E9
b = 2.56E-10
gamma = 0.8

for i in range(len(args.file)):
    f = DADF5(args.file[i],args.geom[i],'a')
    f.addCalculation('tot_density','np.sum(#rho_dip#;#rho_mob#)')
    f.addGradient(['tot_density'])
    f.addNorm(['grad_FFT'])
    f.displacement_data_proc()
    
    #subgrain size calculations
    for inc in f.get_candidates(f.hdf_file,['normFrobenius_grad_FFT']):
        radius = np.ones(len(f.hdf_file[inc + '/' + 'normFrobenius_grad_FFT']))
        for row in range(len(radius)):
            if f.hdf_file[inc + '/' + 'normFrobenius_grad_FFT'][row] > 0.0:
                radius[row] = 2*gamma/(G*(b**2.0)*(f.hdf_file[inc + '/' + 'normFrobenius_grad_FFT'][row]))
            else:
                radius[row] = 1.0 
        f.hdf_file[inc + '/' + 'r_s'] = radius
           
    #cell/material point size
    f.geom_info()
    unit = 1E-06   #user needs to give the units currently (m)
    D    = 50E-06  # average grain size in the material [upper limit for possible subgrain size][better would be to get distance to GB for filtering]
    cell_size = min((f.size/f.grid)[0:2])*unit   #this change is for 2D datasets
    print('cell size is',cell_size)
    neighbour_radius = max((f.size/f.grid)[0:2])*unit
    #ISSUE: if we have a very coarse grid, then we will generate excess nuclei
    # this is because the critical subgrain sizes will always be smaller than the cell size
    # this gives the upper limit to the simulation resolution
    # cell_size more than 50 um is not desired
    # need to check this approach with manual thresholds
    
    #instantaneous nuclei
    for count,inc in enumerate(f.get_candidates(f.hdf_file,['r_s'])):
        rotation_path = f.get_candidates(f.hdf_file,['grain_rotation']) #'/inc00036/constituent/1_low_alloy_steel/grain_rotation'
        f.hdf_file[inc + '/nucleation_flag'] = np.zeros(len(f.hdf_file[inc + '/r_s']))
        for row in range(len(f.hdf_file[inc + '/r_s'])):
            if 2.0*(f.hdf_file[inc + '/r_s'][row]) < cell_size and f.hdf_file[rotation_path[count] + '/grain_rotation'][row][3] > 5.0:
                #change dislocation density of these cells, the orientation is inherited
#                print('row',row)
#                print('r_s',f.hdf_file[inc + '/r_s'][row])
                f.hdf_file[inc + '/nucleation_flag'][row] = 1
                if args.restart:
                    restart_file['PlasticPhases/ 1 _convergedStateConst'][row,0:48] = 1E12 #the number 48 will change with lattice structure,BCC - 48, FCC - 24
               # else:
               #     f.hdf_file[inc + '/rho_dip'][row] = 1E12
               #     f.hdf_file[inc + '/rho_mob'][row] = 1E12

           # elif 2.0*(f.hdf_file[inc + '/r_s'][row]) > cell_size and f.hdf_file[rotation_path[count] + '/grain_rotation'][row][3] > 5.0: 
#          #      print('row_cell',row)
#          #      print('r_s_cell',f.hdf_file[inc + '/r_s'][row])

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# this part is preferable when you are doing phase field
# the idea is okay but needs to be optimised, only the immediate neighbours should be checked
# currently gives problems with elements at the edges, because there, the neighbourhood changes
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        for row in range(len(f.hdf_file[inc + '/r_s'])):
            if f.hdf_file[inc + '/nucleation_flag'][row] == 1:
                pos_inc_list = f.get_candidates(f.hdf_file,['orig_pos'])
                tree = spatial.KDTree(np.array(f.hdf_file[pos_inc_list[count] + '/orig_pos'])*unit)   #tree consisting of positions multiplied by unit
                neighbour_count = tree.query(np.array(f.hdf_file[pos_inc_list[count] + '/orig_pos'])[row]*unit,
                                                        9)[1]
                print(count,row,neighbour_count)
                good_neighbours = 0
                # current_energy = 0.0
                # new_energy     = 0.0
                for neighbour in neighbour_count:
                    if f.hdf_file[inc + '/nucleation_flag'][neighbour] == 1 and neighbour != row:
                        good_neighbours = good_neighbours + 1
                        print(neighbour)
                if good_neighbours > len(neighbour_count)/2:
                    for neighbour in neighbour_count:
                        f.hdf_file[inc + '/nucleation_flag'][neighbour] = 1 
                



#or modifying the restart file 

#check for GB properties. Need to get the GB energy,mobility vs misorientation plot for austenite
                

#modify the geometry file

