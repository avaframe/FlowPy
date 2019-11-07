#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:15:37 2019

@author: neuhauser
This is core function
"""

import numpy as np
import time
from gravi_class import Cell
import sys
import multiprocessing as mp
import psutil

from PyQt5.QtCore import QThread, pyqtSignal


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

def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]
        
def split_release(release, header_release):
    nodata = header_release["noDataValue"]
    release[release == nodata] = 0
    release[release > 1] = 1
    summ = np.sum(release)
    sum_per_split = summ/mp.cpu_count()
    release_list = []
    breakpoint_x = 0
    
    
    for i in range(breakpoint_x, release.shape[1]):
        if len(release_list) == (psutil.cpu_count() -1):
            c = np.zeros_like(release)
            c[:, breakpoint_x:] = release[:,breakpoint_x:]
            release_list.append(c)
            break
        if np.sum(release[:,breakpoint_x:i]) < sum_per_split:
            continue
        else:
            c = np.zeros_like(release)
            c[:, breakpoint_x:i] = release[:,breakpoint_x:i]
            release_list.append(c)
            breakpoint_x = i
        
    return release_list

    
def calculation(args):
    
    dem = args[0]
    header = args[1]
    forest = args[2]
    process = args[3]
# =============================================================================
#     row_list = args[4]
#     col_list = args[5]
# =============================================================================
    release = args[4]
    
    elh = np.zeros_like(dem)
    mass_array = np.zeros_like(dem)
    count_array = np.zeros_like(dem)

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
            row, col, mass, kin_e = cells.calc_distribution()
            if len(mass) > 0:
                # mass, row, col  = list(zip(*sorted(zip( mass, row, col), reverse=False)))
                kin_e, mass, row, col = list(zip(*sorted(zip(kin_e, mass, row, col), reverse=False)))
                # Sort this lists by mass to start the spreading from the middle

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

        release[elh > 0] = 0  # Check if i hited a release Cell, if so set it to zero and get again the indexes of release cells
        # ToDo: if i hit a startcell add this "mass"
        # ToDo: Backcalulation
        row_list, col_list = get_start_idx(dem, release)
        startcell_idx += 1
    end = time.time()            
    print('\n Time needed: ' + str(end - start) + ' seconds')
    #self.quit()
    return elh, mass_array, count_array


class Simulation(QThread):
    value_changed = pyqtSignal(float)
    finished = pyqtSignal(list, list, list)

    def __init__(self, dem, header, release, release_header, forest, process):
        QThread.__init__(self)
        self.dem = dem
        self.header = header
        self.release = release
        self.release_header = release_header
        self.forest = forest
        self.process = process
        self.numberofprocesses = mp.cpu_count()


    def run(self):

        # This part is for Calculation of all release cells
# =============================================================================
#         row_list, col_list = get_start_idx(self.dem, self.release)
#         divided_rowlist = list(divide_chunks(row_list, int(len(row_list)/self.numberofprocesses - 1)))
#         divided_collist = list(divide_chunks(col_list, int(len(col_list)/self.numberofprocesses - 1)))
# 
#         iterable = []
#         for i in range(self.numberofprocesses):
#             iterable.append((self.dem, self.header, self.forest, self.process, divided_rowlist[i], divided_collist[i]))
# =============================================================================
        
        # This part will is for Calculation of the top release cells and ereasing the lower ones
        
        release_list = split_release(self.release, self.release_header)
        iterable = []
        for i in range(self.numberofprocesses):
            iterable.append((self.dem, self.header, self.forest, self.process, release_list[i]))
        
        pool = mp.Pool(processes = self.numberofprocesses)
        results = pool.map(calculation, iterable)
        pool.close()
        pool.join()
        
        print("Processes finished")
            
        elh_list = []
        mass_list = []
        cc_list = []
        for i in range(len(results)):
            res = results[i]
            res = list(res)
            elh_list.append(res[0])
            mass_list.append(res[1])
            cc_list.append(res[2])        

        self.finished.emit(elh_list, mass_list, cc_list)
        print("Results passed")       
        