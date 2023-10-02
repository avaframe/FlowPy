#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:15:37 2019

This is the core function for Flow-Py, it handles: 
- Sorting release pixels by altitude(get_start_idx)
- Splitting function of the release layer for multiprocessing(split_release)
- Back calculation if infrastructure is hit
- Calculation of run out, etc. (Creating the cell_list and iterating through
the release pixels, erasing release pixels that were hit, stop at the border 
of DEM, return arrays)


    Copyright (C) <2020>  <Michael Neuhauser>
    Michael.Neuhauser@bfw.gv.at

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from flow_class import Cell


def get_start_idx(dem, release):
    """Sort Release Pixels by altitude and return the result as lists for the 
    Rows and Columns, starting with the highest altitude
    
    Input parameters:
        dem         Digital Elevation Model to gain information about altitude
        release     The release layer, release pixels need int value > 0
        
    Output parameters:
        row_list    Row index of release pixels sorted by altitude
        col_list    Column index of release pixels sorted by altitude
        """
    row_list, col_list = np.where(release > 0)  # Gives back the indices of the release areas
    if len(row_list) > 0:
        altitude_list = []
        for i in range(len(row_list)):
            altitude_list.append(dem[row_list[i], col_list[i]])    
        altitude_list, row_list, col_list = list(zip(*sorted(zip(altitude_list, row_list, col_list), reverse=True)))
        # Sort this lists by altitude
    return row_list, col_list   


def back_calculation(back_cell):
    """Here the back calculation from a run out pixel that hits a infrastructure
    to the release pixel is performed.
    
    Input parameters:
        hit_cell_list        All cells that hit a Infrastructure
        
    Output parameters:
        Back_list   List of pixels that are on the way to the start cell
                    Maybe change it to array like DEM?
    """
    #start = time.time()
    #if len(hit_cell_list) > 1:
        #hit_cell_list.sort(key=lambda cell: cell.altitude, reverse=False)
        #print("{} Elements sorted!".format(len(hit_cell_list)))
    back_list = []
    for parent in back_cell.parent:
        if parent not in back_list:
            back_list.append(parent)
    for cell in back_list:
        for parent in cell.parent:
            # Check if parent already in list
            if parent not in back_list:
                back_list.append(parent)
    #end = time.time()            
    #print('\n Backcalculation needed: ' + str(end - start) + ' seconds')
    return back_list


def divide_chunks(l, n):
    """Splitting release list in equivalent sub lists, was done before 
    split_release, maybe don't needed anymore... """
    for i in range(0, len(l), n):
        yield l[i:i+n]


def split_release(release, header_release, pieces):
    """Split the release layer in several tiles, the number is depending on the
    available CPU Cores, so every Core gets one tile. The area is determined by
    the number of release pixels in it, so that every tile has the same amount
    of release pixels in it. Splitting in x(Columns) direction. 
    The release tiles have still the size of the original layer, so no split
    for the DEM is needed.
    
    Input parameters: 
        release         the release layer with release pixels as int > 0
        header_release  the header of the release layer to identify the 
                        noDataValue
                        
    Output parameters:
        release_list    A list with the tiles(arrays) in it [array0, array1, ..]
        """
        
    nodata = header_release["noDataValue"]
    if nodata:
        release[release == nodata] = 0
    else:
        print("Release Layer has no No Data Value, negative Value asumed!")
        release[release < 0] = 0
    release[release > 1] = 1
    summ = np.sum(release) # Count number of release pixels
    print("Number of release pixels: ", summ)
    sum_per_split = summ/pieces  # Divide the number by avaiable Cores
    release_list = []
    breakpoint_x = 0

    for i in range(breakpoint_x, release.shape[1]):
        if len(release_list) == (pieces - 1):
            c = np.zeros_like(release)
            c[:, breakpoint_x:] = release[:, breakpoint_x:]
            release_list.append(c)
            break
        if np.sum(release[:, breakpoint_x:i]) < sum_per_split:
            continue
        else:
            c = np.zeros_like(release)
            c[:, breakpoint_x:i] = release[:, breakpoint_x:i]
            release_list.append(c)
            print("Release Split from {} to {}".format(breakpoint_x, i))
            breakpoint_x = i

    return release_list
    
def get_idx_max(energy):   
    #new PAULA
    row_max,col_max = np.where(energy == energy.max())
    row_max = round(np.mean(row_max))
    col_max = round(np.mean(col_max))    
    return (col_max, row_max)

def bresenham_line(start, end):
    #new PAULA
    """
    Generate coordinates of the line between start and end using Bresenham's algorithm.
    """
    x0, y0 = start
    x1, y1 = end
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx - dy
    
    line_coords = []
    
    while True:
        line_coords.append((x0, y0))
        if x0 == x1 and y0 == y1:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x0 += sx
        if e2 < dx:
            err += dx
            y0 += sy
                
    return line_coords    
    
def get_2d_arrays(start, end, raster_3d):
    #new PAULA
    s = bresenham_line(start, end)
    array_2d = [raster_3d[y, x] for x, y in s]
    return s, array_2d


def get_coord_array(dem):
    # NEW PAULA
    num_rows, num_cols = dem.shape
    # Initialize x and y as 2D arrays filled with zeros
    row_idx = [[0 for _ in range(num_cols)] for _ in range(num_rows)]
    col_idx = [[0 for _ in range(num_cols)] for _ in range(num_rows)]

    for i in range(num_rows):
        for j in range(num_cols):
            row_idx[i][j] = i  # Row index
            col_idx[i][j] = j  # Column index

    return row_idx,col_idx
    

def plot_path(startcell_idx, row_idx_start, col_idx_start, dem, flux_array, z_delta_array, energy):
    #new PAULA

    fig, axs = plt.subplots(2,2) 

    fig.set_figheight(10)
    fig.tight_layout(pad=3.0)
    fig.set_figwidth(20)
    
    axs[0,0].imshow(dem, cmap ='Greys', alpha=0.8)
    f = axs[0,0].imshow(energy,  cmap = 'Blues', alpha = 0.6)
    fig.colorbar(f, ax = axs[0,0])
    axs[0,0].set(ylabel = 'y in [m]')
    
    axs[0,1].plot(dem[row_idx,:], label = 'topography, center line')
    axs[0,1].plot(z_delta_array[row_idx,:] + dem[row_idx,:], label = f'Z_delta')
    #axs[0,1].axhline(0)
    axs[0,1].set(ylabel='altitude Z')
    axs[0,1].legend()
    
    axs[1,1].plot(energy[row_idx,:], label='flow energy')
    axs[1,1].set(ylabel='flow energy = flux * z_delta / 2')
    
    x = np.arange(len(dem[row_idx,:]))
    axs[1,0].bar(x,flux_array[row_idx,:])
    axs[1,0].set(ylabel='flux')
    
    i,j = axs.shape
    for idx in np.arange(i):
    	for jdx in np.arange(j):
        	axs[idx,jdx].set(xlabel='x in [m]')
   
        	
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/Z_delta_path_{startcell_idx}.png')

def plot_path_s(startcell_idx, row_idx_start, col_idx_start, cellsize, dem, flux, z_delta, energy):
    '''
    Plot z_delta, flow_energy, flux with x-axis = s(x,y)
    '''
    #new PAULA
    start = (col_idx_start, row_idx_start)
    end = get_idx_max(energy)
    s, dem_2d = get_2d_arrays(start, end, dem)
    s, flux_2d = get_2d_arrays(start, end, flux)
    s, z_delta_2d = get_2d_arrays(start, end, z_delta)
    s, energy_2d = get_2d_arrays(start, end, energy)
    fig, axs = plt.subplots(2,2) 

    fig.set_figheight(10)
    fig.tight_layout(pad=3.0)
    fig.set_figwidth(20)
    
    axs[0,0].imshow(dem, cmap ='Greys', alpha=0.8)
    f = axs[0,0].imshow(energy,  cmap = 'Blues', alpha = 0.6)
    fig.colorbar(f, ax = axs[0,0])
    axs[0,0].set(ylabel = 'y in [m]')
    axs[0,0].legend()
    
    axs[0,1].plot(dem_2d, label = 'topography, center line')
    axs[0,1].plot([d + z for d,z in zip(dem_2d, z_delta_2d)], label = f'Z_delta')
    #axs[0,1].axhline(0)
    axs[0,1].set(ylabel='altitude Z')
    axs[0,1].legend()
    
    axs[1,1].plot(energy_2d, label='flow energy')
    axs[1,1].set(ylabel='flow energy = flux * z_delta / 2')
    
    x = np.arange(len(dem_2d))
    axs[1,0].bar(x,flux_2d)
    axs[1,0].set(ylabel='flux')
    

    i,j = axs.shape
    for idx in np.arange(i):
    	for jdx in np.arange(j):    
	        axs[idx,jdx].set(xlabel='iteration steps')

    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/PATH_S_{startcell_idx}.png')



def calc_thalweg_centerof(generation, co_parameter_path, parameter_path):
    # NEW PAULA    

    parameter_it_sum = np.zeros(int(np.max(generation))+1)
    co_it_sum = np.zeros(int(np.max(generation))+1)
    parameter_coF = np.zeros(int(np.max(generation))+1)
    parameter = np.array(parameter_path)
    co_param = np.array(co_parameter_path)

    for i in np.arange(int(np.max(generation)) + 1):
        idx = np.where(generation == i)

        # sum of parameter for every generation
        parameter_it_sum[i] = np.sum(parameter[generation == i])
        #print(parameter_it_sum[i])
        co_it_sum[i] = np.sum(co_param[idx])

        parameter_coF[i] = 1 / co_it_sum[i] * np.sum(parameter[idx] * co_param[idx])
    '''
    for i in range(10):
        print(i, generation[generation == i])
    

    fig,ax = plt.subplots(2)
    fig.tight_layout(pad=3.0)

    ax[0].plot(parameter_it_sum)
    ax[0].set_ylim(0,1.2)
    ax[0].set_xlabel('generation')
    ax[0].set_ylabel('sum of flux')

    f = ax[1].imshow(generation)
    fig.colorbar(f, ax = ax[1], label = 'generation')

    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/generation2.png')
    '''
    return parameter_it_sum, parameter_coF


def plot_path_sept23(dem,cellsize, generation_path, row, col, flux_path, energy, altitude, travel_length, z_delta):

    #PAULA
    n_startcell = energy.shape[0]

    for n in range(n_startcell): 
        fig, axs = plt.subplots(4,2) 

        fig.set_figheight(10)
        fig.tight_layout(pad=3.0)
        fig.set_figwidth(20)

        variables = {'col': col, 'row': row, 'flux':flux_path[n,:,:], 'energy': energy[n,:,:],
         'altitude': altitude[n,:,:], 's' : travel_length[n,:,:], 'z_delta': z_delta[n,:,:]}
        for var_name, var in variables.items():
            sumF, coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], var) # center of flux of every variable
            sumE, coE = calc_thalweg_centerof(generation_path[n,:,:], energy[n,:,:], var) # center of energy of every variable
    
            globals()[f'{var_name}_sumF'] = sumF
            globals()[f'{var_name}_coF'] = coF
            globals()[f'{var_name}_sumE'] = sumE
            globals()[f'{var_name}_coE'] = coE

        '''
        #center of flux
        row_sumF, row_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], row) # of row position
        col_sumF, col_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], col) # of col position
        flux_sumF, flux_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], flux_path[n,:,:]) # of flux
        altitude_sumF, altitude_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], altitude[n,:,:]) # of altitude difference to startcell
        s_sumF, s_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], travel_length[n,:,:]) # of travel_length
        z_delta_sumF, z_delta_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], z_delta[n,:,:]) # of z_delta
        energy_sumF, energy_coF = calc_thalweg_centerof(generation_path[n,:,:], flux_path[n,:,:], energy[n,:,:]) # of z_delta
        '''

        axs[0,0].imshow(dem, cmap ='Greys', alpha=0.8)
        axs[0,0].contour(dem, levels = 10, colors ='k',linewidths=0.5)
        axs[0,0].scatter(col_coF[::10], row_coF[::10], c = 'k', s = 0.4, label = 'center of flux')
        axs[0,0].scatter(col_coE[::10], row_coE[::10], c = 'r', s = 0.4, label = 'center of energy')
        f = axs[0,0].imshow(energy[n,:,:], cmap = 'Blues', alpha = 0.6)
        fig.colorbar(f, ax = axs[0,0], label = 'flow energy')
        axs[0,0].set(ylabel = 'y in [m]')
        axs[0,0].set(xlabel = 'x in [m]')
        #x_ticks = np.linspace(0, len(dem[0]), 6)  # Adjust the step size as needed
        x_ticks = np.arange(0, len(dem[0]),200) # funktioniert bei datensatz mit 1000 gitterpunkten
        x_tick_labels = [str(round(label * cellsize)) for label in x_ticks] 
        axs[0,0].set_xticks(x_ticks)
        axs[0,0].set_xticklabels(x_tick_labels)   
        #y_ticks = np.linspace(0, len(dem[:,0]), 5)  # Adjust the step size as needed
        y_ticks = np.arange(0, len(dem[:,0]),100) # funktioniert bei datensatz mit 400 gitterpunkten
        y_tick_labels = [str(round(label * cellsize)) for label in y_ticks] 
        axs[0,0].set_yticks(y_ticks)
        axs[0,0].set_yticklabels(y_tick_labels) 
        axs[0,0].legend()

        axs[1,0].plot(s_coF, altitude_coF, 'b--', label = 'topography coF')
        axs[1,0].plot(s_coF,[d + z for d,z in zip(altitude_coF, z_delta_coF)], 'b',label = 'z_delta coF')
        axs[1,0].plot(s_coE, altitude_coE, 'r--', label = 'topography coE')
        axs[1,0].plot(s_coE, [d + z for d,z in zip(altitude_coE, z_delta_coE)], 'r',label = 'z_delta coE')
        axs[1,0].set(xlabel = 's in [m]')      
        axs[1,0].set(ylabel = 'altitude in [m]')
        axs[1,0].legend()

        f = axs[2,0].imshow(z_delta[n,:,:], cmap = 'Blues')
        fig.colorbar(f, ax = axs[2,0], label = 'z_delta')

        axs[0,1].plot(altitude_coF, 'b--', label = 'topography coF')
        axs[0,1].plot(altitude_coF + z_delta_coF, 'b',label = 'z_delta coF')
        axs[0,1].plot(altitude_coE, 'r--', label = 'topography coE')
        axs[0,1].plot(altitude_coE + z_delta_coE, 'r', label = 'z_delta coE')
        axs[0,1].set(ylabel = 'altitude Z (coF)')
        axs[0,1].set(xlabel = 'iteration step')  
        axs[0,1].legend()
 
        axs[1,1].plot(flux_sumF, 'b', label = 'coF')
        axs[1,1].plot(flux_sumE, 'r', label = 'coE')
        axs[1,1].set(ylabel = 'sum of flux')
        axs[1,1].set_ylim(0,1.1) 

        axs[2,1].plot(energy_sumF, 'b', label = 'coF')
        axs[2,1].plot(energy_sumE, 'r', label = 'coE')
        axs[2,1].set(ylabel = 'sum of flow energy')

        axs[3,1].plot(s_coF, 'b', label = 'coF')
        axs[3,1].plot(s_coE, 'r', label = 'coE')
        axs[3,1].set(ylabel = 'travel length')

        #for all axes on right side
        for i in range(4):
            axs[i,1].legend()
            axs[i,1].set(xlabel = 'iteration step')   
        
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/path_analysis_startcell{n}.png')

    	
    
def calculation(args):
    """This is the core function where all the data handling and calculation is
    done. 
    
    Input parameters:
        dem         The digital elevation model
        header      The header of the elevation model
        infras      The infrastructure layer
        process     Which process to calculate (Avalanche, Rockfall, SoilSlides)     
        release     The list of release arrays
        alpha
        exp
        flux_threshold
        max_z_delta
        
    Output parameters:
        z_delta_array   Array like DEM with the max. Energy Line Height for every 
                        pixel
        z_delta_sum     Array...
        mass_array  Array with max. concentration factor saved
        count_array Array with the number of hits for every pixel
        elh_sum     Array with the sum of Energy Line Height
        back_calc   Array with back calculation, still to do!!!
        """
    
    dem = args[0]
    header = args[1]
    infra = args[2]
    forest = args[3]
    release = args[4]
    alpha = args[5]
    exp = args[6]
    flux_threshold = args[7]
    max_z_delta = args[8]
    #print(len(args), max_z_delta)
    
    # JT PROGRAMMIERT CRAZY STUFF
    #ch alti
    altitude_diff_array = np.zeros_like(dem)
    #ch alti
    travel_length_array = np.zeros_like(dem)
    # JT ENDE
    z_delta_array = np.zeros_like(dem)
    z_delta_sum = np.zeros_like(dem)
    flux_array = np.zeros_like(dem)
    flow_energy_array = np.zeros_like(dem)
    count_array = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)
    fp_travelangle_array = np.zeros_like(dem)
    sl_travelangle_array = np.ones_like(dem) * 90
    back_list = []

    cellsize = header["cellsize"]
    nodata = header["noDataValue"]

    # Core
    start = datetime.now().replace(microsecond=0)
    row_list, col_list = get_start_idx(dem, release)

    startcell_idx = 0
    while startcell_idx < len(row_list):
        
        sys.stdout.write('\r' "Calculating Startcell: " + str(startcell_idx + 1) + " of " + str(len(row_list)) + " = " + str(
            round((startcell_idx + 1) / len(row_list) * 100, 2)) + "%" '\r')
        sys.stdout.flush()

        cell_list = []
        row_idx = row_list[startcell_idx]
        col_idx = col_list[startcell_idx]
        dem_ng = dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2]  # neighbourhood DEM
        if (nodata in dem_ng) or np.size(dem_ng) < 9:
            startcell_idx += 1
            continue

        startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1, 0, None,
                         alpha, exp, flux_threshold, max_z_delta, startcell=True)
        # If this is a startcell just give a Bool to startcell otherwise the object startcell

        cell_list.append(startcell)

        for idx, cell in enumerate(cell_list):
            row, col, flux, z_delta = cell.distribution()

            if len(flux) > 0:
                # mass, row, col  = list(zip(*sorted(zip( mass, row, col), reverse=False)))
                
                z_delta, flux, row, col = list(zip(*sorted(zip(z_delta, flux, row, col), reverse=False)))
                # Sort this lists by elh, to start with the highest cell

            for i in range(idx, len(cell_list)):  # Check if Cell already exists
                k = 0
                while k < len(row):
                    if row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex:
                        cell_list[i].add_os(flux[k])
                        cell_list[i].add_parent(cell)
                        if z_delta[k] > cell_list[i].z_delta:
                            cell_list[i].z_delta = z_delta[k]
                        row = np.delete(row, k)
                        col = np.delete(col, k)
                        flux = np.delete(flux, k)
                        z_delta = np.delete(z_delta, k)
                    else:
                        k += 1

            for k in range(len(row)):
                dem_ng = dem[row[k] - 1:row[k] + 2, col[k] - 1:col[k] + 2]  # neighbourhood DEM
                if (nodata in dem_ng) or np.size(dem_ng) < 9:
                    continue
                cell_list.append(
                    Cell(row[k], col[k], dem_ng, cellsize, flux[k], z_delta[k], cell, alpha, exp, flux_threshold, max_z_delta, startcell))
	    # JT PROGRAMMIERT CRAZY STUFF
        # ch alti
            altitude_diff_array[cell.rowindex, cell.colindex] = max(altitude_diff_array[cell.rowindex, cell.colindex], cell.altitude_diff) 
        #ch alti
            travel_length_array[cell.rowindex, cell.colindex] = max(travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
	    # JT ENDE
            z_delta_array[cell.rowindex, cell.colindex] = max(z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
            flux_array[cell.rowindex, cell.colindex] = max(flux_array[cell.rowindex, cell.colindex], cell.flux)
            count_array[cell.rowindex, cell.colindex] += 1
            z_delta_sum[cell.rowindex, cell.colindex] += cell.z_delta
            fp_travelangle_array[cell.rowindex, cell.colindex] = max(fp_travelangle_array[cell.rowindex, cell.colindex], cell.max_gamma)
            sl_travelangle_array[cell.rowindex, cell.colindex] = max(sl_travelangle_array[cell.rowindex, cell.colindex], cell.sl_gamma)
            #PAULA
            flow_energy_array[cell.rowindex, cell.colindex] = max(flow_energy_array[cell.rowindex, cell.colindex], cell.flow_energy)
            #ende paula
            
        #Backcalculation
            if infra[cell.rowindex, cell.colindex] > 0:
                #backlist = []
                back_list = back_calculation(cell)

                for back_cell in back_list:
                    backcalc[back_cell.rowindex, back_cell.colindex] = max(backcalc[back_cell.rowindex, back_cell.colindex],
                                                                           infra[cell.rowindex, cell.colindex])
        release[z_delta_array > 0] = 0
        # Check if i hit a release Cell, if so set it to zero and get again the indexes of release cells
        row_list, col_list = get_start_idx(dem, release)
        startcell_idx += 1
    end = datetime.now().replace(microsecond=0)
    #elh_multi[elh_multi == 1] = 0         
    print('\n Time needed: ' + str(end - start))
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array
    #JT CRAZY STUFF
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array
    #JT ENDE
    #ch alti
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array, altitude_diff_array
    #ch alti
    #PAULA
    return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array, altitude_diff_array,flow_energy_array
    #ende Paula
    
def calculation_effect(args):
    """This is the core function where all the data handling and calculation is
    done. 
    
    Input parameters:
        dem         The digital elevation model
        header      The header of the elevation model
        process     Which process to calculate (Avalanche, Rockfall, SoilSlides)     
        release     The list of release arrays
        
    Output parameters:
        elh         Array like DEM with the max. Energy Line Height for every 
                    pixel
        mass_array  Array with max. concentration factor saved
        count_array Array with the number of hits for every pixel
        elh_sum     Array with the sum of Energy Line Height
        back_calc   Array with back calculation, still to do!!!
        """
    
    dem = args[0]
    header = args[1]
    forest = args[2]
    release = args[3]
    alpha = args[4]
    exp = args[5]
    flux_threshold = args[6]
    max_z_delta = args[7]

    # CH PROGRAMMIERT CRAZY STUFF
    #ch alti
    altitude_diff_array = np.zeros_like(dem)
    #ch alti
    travel_length_array = np.zeros_like(dem)
    # CH ENDE
    z_delta_array = np.zeros_like(dem)
    z_delta_sum = np.zeros_like(dem)
    flux_array = np.zeros_like(dem)
    flow_energy_array = np.zeros_like(dem)
    count_array = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)
    fp_travelangle_array = np.zeros_like(dem)  # fp = Flow Path
    sl_travelangle_array = np.zeros_like(dem)  # sl = Straight Line

    cellsize = header["cellsize"]
    nodata = header["noDataValue"]

    # Core
    start = datetime.now().replace(microsecond=0)
    row_list, col_list = get_start_idx(dem, release)
    
    #PAULA
    flux_path = np.zeros((len(row_list), dem.shape[0],dem.shape[1]))
    flow_energy_path = np.zeros((len(row_list), dem.shape[0],dem.shape[1]))
    z_delta_path = np.zeros_like(flux_path)
    travel_length_path = np.zeros_like(flux_path)
    generation_path = np.zeros_like(flux_path)
    altitude_path = np.zeros_like(flux_path)
    #ende paula
    	
    startcell_idx = 0
    while startcell_idx < len(row_list):
        #calculate for every startcell
        
        sys.stdout.write('\r' "Calculating Startcell: " + str(startcell_idx + 1) + " of " + str(len(row_list)) + " = " + str(
            round((startcell_idx + 1) / len(row_list) * 100, 2)) + "%" '\r')
        sys.stdout.flush()

        cell_list = []
        row_idx = row_list[startcell_idx]
        col_idx = col_list[startcell_idx]
        
        dem_ng = dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2]  # neighbourhood DEM
        if (nodata in dem_ng) or np.size(dem_ng) < 9:
            startcell_idx += 1
            continue
        startcell = Cell(row_idx, col_idx, dem_ng, forest[row_idx, col_idx], 
                         cellsize, 1, 0, None,
                         alpha, exp, flux_threshold, max_z_delta, True)
        
        # If this is a startcell just give a Bool to startcell otherwise the object startcell

        cell_list.append(startcell)
        
    
            
        for idx, cell in enumerate(cell_list):
          
            row, col, flux, z_delta = cell.calc_distribution()

            if len(flux) > 0:
                z_delta, flux, row, col = list(zip(*sorted(zip(z_delta, flux, row, col), reverse=False)))  # reverse = True == descending
            
            for i in range(idx, len(cell_list)):  # Check if Cell already exists
                k = 0
                while k < len(row):
                    if row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex:
                        cell_list[i].add_os(flux[k])
                        cell_list[i].add_parent(cell)
                        if z_delta[k] > cell_list[i].z_delta:
                            cell_list[i].z_delta = z_delta[k]

                        row = np.delete(row, k)
                        col = np.delete(col, k)
                        flux = np.delete(flux, k)
                        z_delta = np.delete(z_delta, k)
                    else:
                        k += 1
            
          
                            

            for k in range(len(row)):
                dem_ng = dem[row[k] - 1:row[k] + 2, col[k] - 1:col[k] + 2]  # neighbourhood DEM
                if (nodata in dem_ng) or np.size(dem_ng) < 9:
                    continue
                cell_list.append(
                     Cell(row[k], col[k], dem_ng, forest[row[k], col[k]],
                          cellsize, flux[k], z_delta[k], cell, alpha, exp, 
                          flux_threshold, max_z_delta, startcell))
                
                



        
        for cell in cell_list:
        # Write values for every cell in an array 
        # For raster raster level: Overlay of all path-results (of all startcells) (2 dim)
        # For path level: Write path results (for each startcells) in a dimension of the array (3 dim)
        #PATH EBENE
        # ch PROGRAMMIERT CRAZY STUFF
            #ch alti
            #altitude_diff_array[cell.rowindex, cell.colindex] = max(altitude_diff_array[cell.rowindex, cell.colindex], cell.startcell.altitude - cell.altitude) 
            altitude_diff_array[cell.rowindex, cell.colindex] = max(altitude_diff_array[cell.rowindex, cell.colindex], cell.altitutde_diff) 
            #ch alti
            travel_length_array[cell.rowindex, cell.colindex] = max(travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
        # CH ENDE
            z_delta_array[cell.rowindex, cell.colindex] = max(z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
            
            #PAULA
            flow_energy_array[cell.rowindex, cell.colindex] = max(flow_energy_array[cell.rowindex, cell.colindex], cell.flow_energy)

            # arrays for path 
            flux_path[startcell_idx,cell.rowindex, cell.colindex] = max(flux_path[startcell_idx,cell.rowindex, cell.colindex],cell.flux)
            flow_energy_path[startcell_idx,cell.rowindex, cell.colindex] = max(flux_path[startcell_idx,cell.rowindex, cell.colindex],cell.flow_energy)
            z_delta_path[startcell_idx,cell.rowindex, cell.colindex] = max(z_delta_path[startcell_idx,cell.rowindex, cell.colindex],cell.z_delta)
            travel_length_path[startcell_idx,cell.rowindex, cell.colindex] = max(travel_length_path[startcell_idx,cell.rowindex, cell.colindex],cell.min_distance)
            generation_path[startcell_idx,cell.rowindex, cell.colindex] = max(generation_path[startcell_idx,cell.rowindex, cell.colindex],cell.generation)
            altitude_path[startcell_idx,cell.rowindex, cell.colindex] = max(altitude_path[startcell_idx,cell.rowindex, cell.colindex], cell.altitude) 

            #ENDE Paula
            
            flux_array[cell.rowindex, cell.colindex] = max(flux_array[cell.rowindex, cell.colindex],
                                                        cell.flux)
            count_array[cell.rowindex, cell.colindex] += 1
            z_delta_sum[cell.rowindex, cell.colindex] += cell.z_delta
            fp_travelangle_array[cell.rowindex, cell.colindex] = max(fp_travelangle_array[cell.rowindex, cell.colindex],
                                                                    cell.max_gamma)
            sl_travelangle_array[cell.rowindex, cell.colindex] = max(sl_travelangle_array[cell.rowindex, cell.colindex],
                                                                    cell.sl_gamma)
        
                                                                    
	#PAULA
        #plot_path(startcell_idx, row_idx, col_idx, dem, flux_array, z_delta_array, flow_energy_array)
        plot_path_s(startcell_idx, row_idx, col_idx, cellsize, dem, flux_array, z_delta_array, flow_energy_array)   
        #ende paula
        startcell_idx += 1

        


    #PAULA
    row_idx, col_idx = get_coord_array(dem)
    plot_path_sept23(dem, cellsize, generation_path, row_idx, col_idx, flux_path, flow_energy_path, altitude_path, travel_length_path, z_delta_path)


    #oaula ende
    end = datetime.now().replace(microsecond=0)        
    print('\n Time needed: ' + str(end - start))
    #CH 
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array
    #CH ENDE
    #ch alti
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array, altitude_diff_array, flux_path
    #ch alti
    #PAULA
    return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array, altitude_diff_array, flux_path, flow_energy_array
    #ende paula
    
