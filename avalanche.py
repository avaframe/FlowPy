# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:44:52 2019

@author: Michael Neuhauser
This is model to simulate gravitave processes.
alpha = ...
exp = 8 for avalanche; 75 for single flow like rockfall
"""

import numpy as np
import math

class cell():
    
    def __init__(self, rowindex, colindex, altitude):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = altitude
        self.dx = None
        self.dy = None
        self.dz = None
        self.cellsize = 10
        self.exp = 8
               
    def calc_vector(self, startcell):
        self.dx = np.abs(startcell.colindex - self.colindex)
        self.dy = np.abs(startcell.rowindex - self.rowindex)
        self.dz = np.abs(startcell.altitude - self.altitude)
        
    def calc_distribution(self, dem_ng):
        self.tan_beta = np.zeros((3, 3))
        self.mass_dist = np.zeros((3,3))
        exp = self.exp
        cellsize = self.cellsize
        threshold = 0.05

        for i in range(-1, 2):
            for j in range(-1, 2):

                row_idx = self.rowindex + i
                col_idx = self.colindex + j
                #height = dem[row_idx, col_idx]
                if row_idx == self.rowindex or col_idx == self.colindex:
                    distance = cellsize
                else:
                    distance = cellsize * math.sqrt(2)
                self.tan_beta[i+1, j+1] = (dem_ng[i+1, j+1] - self.altitude)/distance

        self.tan_beta[self.tan_beta > 0] = 0
        for i in range(3):
            for j in range(3):
                self.mass_dist[i, j] = self.tan_beta[i, j] ** exp / np.sum(self.tan_beta ** exp)
        out_row_local, out_col_local = np.where(self.mass_dist > threshold)
        out_row = self.rowindex - 1 + out_row_local
        out_col = self.colindex - 1 + out_col_local
        return out_row, out_col
    
    
dem = np.ones((10, 10), dtype=np.int)
mass_dist = np.zeros_like(dem)

for i in range(10):
    dem[:,i] = dem[:,i] * (10-i)

cell_list = []  
startcell = cell(5, 1, dem[5, 1])
mass_dist[startcell.rowindex, startcell.colindex] = 1

i, j = startcell.calc_distribution(dem[startcell.rowindex-1:startcell.rowindex+2,startcell.colindex-1:startcell.colindex+2])

for k in range(len(i)):
    cell_list.append(cell(i[k], j[k], dem[i[k], j[k]]))
    mass_dist[i[k], j[k]] = 1
    
for cells in cell_list:    
    i, j = cells.calc_distribution(dem[cells.rowindex-1:cells.rowindex+2,cells.colindex-1:cells.colindex+2])    
    k = 0    
    while k < len(i): # Controll if cell is on boarder of DEM, when delete it from the list
        if i[k] == np.size(dem, 0) - 1 or j[k] == np.size(dem, 1) - 1:
            i = np.delete(i, k)
            j = np.delete(j, k)
        k += 1
    for k in range(len(i)): 
        if len(i) > 0:
            cell_list.append(cell(i[k], j[k], dem[i[k], j[k]]))
            mass_dist[i[k], j[k]] = 1
    
    

