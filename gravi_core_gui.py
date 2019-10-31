#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:15:37 2019

@author: neuhauser
This is core function
"""

import numpy as np
from gravi_class import Cell
import sys
from PyQt5.QtCore import QThread, pyqtSignal
import time


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


class Simulation(QThread):
    value_changed = pyqtSignal(float)
    finished = pyqtSignal(np.ndarray, np.ndarray, np.ndarray)

    def __init__(self, dem, header, release, forest, process):
        QThread.__init__(self)
        self.dem = dem
        self.header = header
        self.release = release
        self.forest = forest
        self.process = process

    def run(self):
        elh = np.zeros_like(self.dem)
        mass_array = np.zeros_like(self.dem)
        count_array = np.zeros_like(self.dem)

        cellsize = self.header["cellsize"]
        nodata = self.header["noDataValue"]

        # Core

        row_list, col_list = get_start_idx(self.dem, self.release)
        # =============================================================================
        # for i in range(len(row_list)):
        #     calc_list.append((row_list[i], col_list[i]))
        #
        #
        # cpu_number = len(os.sched_getaffinity(0)) # counts number of cpu free.
        # pool = Pool(processes=4)
        # pool.map(calculation, calc_list)
        # =============================================================================
        start = time.time()
        startcell_idx = 0
        while startcell_idx < len(row_list):
            
            sys.stdout.write('\r' "Calculating Startcell: " + str(startcell_idx + 1) + " of " + str(len(row_list)) + " = " + str(
                round((startcell_idx + 1) / len(row_list) * 100, 2)) + "%" '\r')
            sys.stdout.flush()

            calculation_percent = round((startcell_idx + 1) / len(row_list) * 100, 2)
            self.value_changed.emit(calculation_percent)

            cell_list = []
            row_idx = row_list[startcell_idx]
            col_idx = col_list[startcell_idx]
            dem_ng = self.dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2]  # neighbourhood DEM
            if (nodata in dem_ng) or np.size(dem_ng) < 9:
                startcell_idx += 1
                continue

            startcell = Cell(self.process, row_idx, col_idx, dem_ng, cellsize, 1, 0, self.forest[row_idx, col_idx], None,
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
                    dem_ng = self.dem[row[k] - 1:row[k] + 2, col[k] - 1:col[k] + 2]  # neighbourhood DEM
                    if (nodata in dem_ng) or np.size(dem_ng) < 9:
                        continue
                    cell_list.append(
                        Cell(self.process, row[k], col[k], dem_ng, cellsize, mass[k], kin_e[k], self.forest[row[k], col[k]],
                             cells, startcell))

                elh[cells.rowindex, cells.colindex] = max(elh[cells.rowindex, cells.colindex], cells.kin_e)
                mass_array[cells.rowindex, cells.colindex] = max(mass_array[cells.rowindex, cells.colindex], cells.mass)
                count_array[cells.rowindex, cells.colindex] += 1

            self.release[elh > 0] = 0  # Check if i hited a release Cell, if so set it to zero and get again the indexes of release cells
            # ToDo: if i hit a startcell add this "mass"
            # ToDo: Backcalulation
            row_list, col_list = get_start_idx(self.dem, self.release)
            startcell_idx += 1
        self.finished.emit(elh, mass_array, count_array)
        end = time.time()            
        print('Time needed: ' + str(end - start) + ' seconds')
        #return elh, mass_array, count_array

