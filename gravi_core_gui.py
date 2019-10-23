#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:15:37 2019

@author: neuhauser
This is core function
"""

import numpy as np
import sys
import time
#from multiprocessing import Pool
#import os
import raster_io as io
from gravi_class import Cell


## Programm:        
def get_start_idx(dem, release):
    row_list, col_list = np.where(release > 0)  # Gives back the indices of the release areas
    if len(row_list) > 0:
        altitude_list = []
        for i in range(len(row_list)):
            altitude_list.append(dem[row_list[i], col_list[i]])    
        altitude_list, row_list, col_list = list(zip(*sorted(zip(altitude_list, row_list, col_list), reverse=True)))  #Sort this lists by altitude
    return row_list, col_list   


def back_calculation(cell):
    back_list = []
    for parent in cell.parent:
        back_list.append(parent)
    for cell in back_list:
        for parent in cell.parent:
            back_list.append(parent)
    return back_list


def calculation(dem, forest, process, cellsize, nodata, row_list, col_list, startcell_idx, elh, mass_array, count_array):


    cell_list = []
    row_idx = row_list[startcell_idx]
    col_idx = col_list[startcell_idx]
    dem_ng = dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2]  # neighbourhood DEM
    if (nodata in dem_ng) or np.size(dem_ng) < 9:
        #startcell_idx += 1
        return

    startcell = Cell(process, row_idx, col_idx, dem_ng, cellsize, 1, 0, forest[row_idx, col_idx], None, startcell=True)
    # If this is a startcell just give a Bool to startcell otherwise the object startcell

    cell_list.append(startcell)
    for cells in cell_list:
        row, col, mass, kin_e = cells.calc_distribution()
        if len(mass) > 0:
            # mass, row, col  = list(zip(*sorted(zip( mass, row, col), reverse=False)))
            kin_e, mass, row, col = list(zip(*sorted(zip(kin_e, mass, row, col), reverse=False)))
            # Sort this lists by elh to start from top to bottom

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
                Cell(process, row[k], col[k], dem_ng, cellsize, mass[k], kin_e[k], forest[row[k], col[k]], cells, startcell))

        elh[cells.rowindex, cells.colindex] = max(elh[cells.rowindex, cells.colindex], cells.kin_e)
        mass_array[cells.rowindex, cells.colindex] = max(mass_array[cells.rowindex, cells.colindex], cells.mass)
        count_array[cells.rowindex, cells.colindex] += 1
        return elh, mass_array, count_array
