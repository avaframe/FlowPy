# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:44:52 2019

@author: Michael
"""

import numpy as np
import math

class cell:
    
    def __init__(self):
        self.rowindex = None
        self.colindex = None
        self.altitude = None
        self.dx = None
        self.dy = None
        self.dz = None
        
    def set_parameters(self, rowindex, colindex, altitude):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = altitude
        
    def calc_vector(self, startcell):
        self.dx = np.abs(startcell.colindex - self.colindex)
        self.dy = np.abs(startcell.rowindex - self.rowindex)
        self.dz = np.abs(startcell.altitude - self.altitude)
        
    def calc_distribution(self, dem_ng):
        tan_beta = np.zeros((3, 3))
        mass_dist = np.zeros((3,3))
        exp = 8
        cellsize = 10

        for i in range(-1, 2):
            for j in range(-1, 2):

                row_idx = self.rowindex + i
                col_idx = self.colindex + j
                #height = dem[row_idx, col_idx]
                if row_idx == self.rowindex or col_idx == self.colindex:
                    distance = cellsize
                else:
                    distance = cellsize * math.sqrt(2)
                tan_beta[i+1, j+1] = (dem_ng[i+1, j+1] - self.altitude)/distance
        print(tan_beta)

        tan_beta[tan_beta > 0] = 0
        for i in range(3):
            for j in range(3):
                mass_dist[i, j] = tan_beta[i, j] ** exp / np.sum(tan_beta ** exp)
        print(mass_dist)
        out_row, out_col = np.where(mass_dist > 0)
        return out_row, out_col
    
    

dem = np.ones((10, 10), dtype=np.int)

for i in range(10):
    dem[:,i] = dem[:,i] * (10-i)

cell_list = []

   
startcell = cell()
startcell.set_parameters(1, 1, dem[1, 1])
#cell_list.append(startcell)
i, j = startcell.calc_distribution(dem[0:3,0:3])

for k in range(len(i)):
    new_cell = cell()
    cell_list.append(new_cell.set_parameters(i[k], j[k], dem[i[k], j[k]]))
    
iteration = 0
#while iteration < 3: 
for cell in cell_list:
    i, j = cell.calc_distribution(dem[0:3,0:3])
    for k in range(len(i)):
        new_cell = cell()
        cell_list.append(new_cell.set_parameters(i[k], j[k], dem[i[k], j[k]]))
    
    

# =============================================================================
# new_cell = cell()
# new_cell.set_parameters(5, 5, dem[5, 5])
# new_cell.calc_vector(startcell)
#     
# cell_list = []
# cell_list.append(cell)
# 
# for i in range(10):
#     new_cell = cell()
#     new_cell.set_parameters(i, i, dem[i, i])
#     new_cell.calc_vector()
#     cell_list.append(new_cell)
# =============================================================================

