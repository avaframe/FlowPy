#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:15:37 2019

@author: Michael Neuhauser
This is the core function for Flow-Py, it handles: 
- Sorting release pixels by altitude(get_start_idx)
- Splitting function of the release layer for multiprocessing(split_release)
- Back calculation if infrastructure is hit (still to check if it works)
- Calculation of run out, etc. (Creating the cell_list and iterating through
the release pixels, erasing release pixels that were hit, stop at the border 
of DEM, return arrays)
"""

import sys
import multiprocessing as mp
import numpy as np
import time
from gravi_class import Cell


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
        altitude_list, row_list, col_list = list(zip(*sorted(zip(altitude_list, row_list, col_list), reverse=True)))  #Sort this lists by altitude
    return row_list, col_list   


def back_calculation(hit_cell_list):
    """Here the back calculation from a run out pixel that hits a infrastructure
    to the release pixel is performed.
    
    Input parameters:
        hit_cell_list        All cells that hit a Infrastructure
        
    Output parameters:
        Back_list   List of pixels that are on the way to the start cell
                    Maybe change it to array like DEM?
    """
    #start = time.time()
    if len(hit_cell_list) > 1:
        hit_cell_list.sort(key=lambda cell: cell.altitude, reverse=False)
        #print("{} Elements sorted!".format(len(hit_cell_list)))
    back_list = []
    for cell in hit_cell_list:
        if cell not in back_list:
            for parent in cell.parent:
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


def split_release(release, header_release):
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
    release[release == nodata] = 0
    release[release > 1] = 1
    summ = np.sum(release) # Count number of release pixels
    sum_per_split = summ/mp.cpu_count()  # Divide the number by avaible Cores
    release_list = []
    breakpoint_x = 0

    for i in range(breakpoint_x, release.shape[1]):
        if len(release_list) == (mp.cpu_count() - 1):
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
            breakpoint_x = i
        
    return release_list

    
def calculation(args):
    """This is the core function where all the data handling and calculation is
    done. 
    
    Input parameters:
        dem         The digital elevation model
        header      The header of the elevation model
        forest      The forest layer
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
    infra = args[2]
    forest = args[3]
    process = args[4]
# =============================================================================
#     row_list = args[4]
#     col_list = args[5]
# =============================================================================
    release = args[5]
    
    elh = np.zeros_like(dem)
    elh_sum = np.zeros_like(dem)
    mass_array = np.zeros_like(dem)
    count_array = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)

    cellsize = header["cellsize"]
    nodata = header["noDataValue"]

    # Core
    start = time.time()
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

        startcell = Cell(process, row_idx, col_idx, dem_ng, cellsize, 1, 0, forest[row_idx, col_idx], None,
                         startcell=True)
        # If this is a startcell just give a Bool to startcell otherwise the object startcell

        cell_list.append(startcell)
        for cells in cell_list:
# =============================================================================
#             if (release[cells.rowindex, cells.colindex] > 0 and cells.is_start != True and process == 'Avalanche'):
#                 cells.add_mass(1)
#                 # Adds the mass of a release cell which is then erased on the way down
#                 # So the mass threshold starts later... 
#                 # Works only for release areas not for release pixels!
# =============================================================================
            # First attempt to include back calculation, second try after all
            # cells are created
# =============================================================================
#             if infra[cells.rowindex, cells.colindex] > 0:
#                 back_list = back_calculation(cells) # would be cool if it gives already array back...
#                 for cell in back_list:
#                     backcalc[cell.rowindex, cell.colindex] = 1
# =============================================================================
                
            row, col, mass, kin_e = cells.calc_distribution()
            if len(mass) > 0:
                # mass, row, col  = list(zip(*sorted(zip( mass, row, col), reverse=False)))
                kin_e, mass, row, col = list(zip(*sorted(zip(kin_e, mass, row, col), reverse=False)))
                # Sort this lists by elh, to start with the highest cell

            for i in range(len(cell_list)):  # Check if Cell already exists
                k = 0
                while k < len(row):
                    if row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex:
                        cell_list[i].add_mass(mass[k])
                        cell_list[i].add_parent(cells)
                        cell_list[i].kin_e = max(cell_list[i].kin_e, kin_e[k])
                        row = np.delete(row, k)
                        col = np.delete(col, k)
                        mass = np.delete(mass, k)
                        kin_e = np.delete(kin_e, k)
                    else:
                        k += 1

            for k in range(len(row)):
                dem_ng = dem[row[k] - 1:row[k] + 2, col[k] - 1:col[k] + 2]  # neighbourhood DEM
                if (nodata in dem_ng) or np.size(dem_ng) < 9:
                    continue
                cell_list.append(
                    Cell(process, row[k], col[k], dem_ng, cellsize, mass[k], kin_e[k], forest[row[k], col[k]],
                         cells, startcell))

            elh[cells.rowindex, cells.colindex] = max(elh[cells.rowindex, cells.colindex], cells.kin_e)
            mass_array[cells.rowindex, cells.colindex] = max(mass_array[cells.rowindex, cells.colindex], cells.mass)
            count_array[cells.rowindex, cells.colindex] += 1
            elh_sum[cells.rowindex, cells.colindex] += cells.kin_e
            
        #Backcalculation
        hit_cell_list = []
        for cell in cell_list:
            if infra[cell.rowindex, cell.colindex] > 0:
                hit_cell_list.append(cell)
        if len(hit_cell_list) > 0:
            back_list = back_calculation(hit_cell_list)
            for cell in back_list:
                backcalc[cell.rowindex, cell.colindex] = 1
            
        release[elh > 0] = 0
        # Check if i hited a release Cell, if so set it to zero and get again the indexes of release cells
        # ToDo: if i hit a startcell add this "mass"
        # ToDo: Backcalculation
        row_list, col_list = get_start_idx(dem, release)
        startcell_idx += 1
    end = time.time()            
    print('\n Time needed: ' + str(end - start) + ' seconds')
    # self.quit()
    return elh, mass_array, count_array, elh_sum, backcalc
