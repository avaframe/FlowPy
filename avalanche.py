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
import rasterio
import sys
import time

sys.path.append('/home/neuhauser/Documents/Milan')
import raster_io as io


class Cell():

    
    def __init__(self, rowindex, colindex, altitude, velocity_sqr, mass):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = altitude
        self.mass = mass
        self.velocity_sqr = velocity_sqr
        self.cellsize = 10
        self.exp = 1
        self.alpha = 5
        self.mu_l = np.tan(np.deg2rad(30))
        #self.mu_l = 0
        self.mu_g = np.tan(np.deg2rad(5))
        self.vel_sqr_out = np.zeros((3,3))
        self.kin = np.zeros((3,3))

        
    def calc_tanbeta(self, dem_ng):       
        self.tan_beta = np.zeros((3, 3))
        cellsize = self.cellsize
        for i in range(-1, 2):
            for j in range(-1, 2):
                row_idx = self.rowindex + i
                col_idx = self.colindex + j
    
                if row_idx == self.rowindex or col_idx == self.colindex:
                    distance = cellsize
                else:
                    distance = cellsize * math.sqrt(2)
                self.tan_beta[i+1, j+1] = (dem_ng[i+1, j+1] - self.altitude)/distance

                
    def calc_mass(self, dem_ng):      
        self.calc_tanbeta(dem_ng)
        self.mass_dist = np.zeros((3, 3))
        exp = self.exp
        mass_threshold = 0.0001

        self.tan_beta[self.vel_sqr_out < 0] = 0
        for i in range(3):
            for j in range(3):
                self.mass_dist[i, j] = self.mass * self.tan_beta[i, j] ** exp / np.sum(self.tan_beta ** exp)
        out_row_local, out_col_local = np.where(self.mass_dist > mass_threshold)
        out_row = self.rowindex - 1 + out_row_local
        out_col = self.colindex - 1 + out_col_local
        mass_list = self.mass_dist[out_row_local, out_col_local]
        
                       
    def calc_velocity(self, startcell, dem_ng):
        g = 9.81
        v_sqr = np.zeros((3, 3))
        
        for i in range(-1, 2):
            for j in range(-1, 2):
                if not (i == 0 and j == 0):
                    dx = np.abs(startcell.colindex - self.colindex + j) * self.cellsize
                    dy = np.abs(startcell.rowindex - self.rowindex + i) * self.cellsize
                    ds_global = np.sqrt(dx**2 + dy**2)
                    #print(dx, dy, ds_global,i,j)
                    dx = j * self.cellsize
                    dy = i * self.cellsize
                    ds_local = np.sqrt(dx**2 + dy**2)
                    #print(dx, dy, ds_local)
                    
                    v_sqr_local = self.velocity_sqr + 2 * g * ((self.altitude - dem_ng[i+1, j+1]) - self.mu_l * ds_local)
                    v_sqr_global = 2 * g*((startcell.altitude - dem_ng[i+1, j+1]) - self.mu_g * ds_global)
                    
                    #print('Local: ' + str(v_sqr_local) + '  Global:' + str(v_sqr_global))
       
                    v_sqr[i+1, j+1] = v_sqr_global
                    self.vel_sqr_out = v_sqr
            
            self.calc_mass(dem_ng)
            self.kin = 0.5 * self.mass_dist * self.vel_sqr_out
       
        #v_sqr[v_sqr < 0] = 0
        return v_sqr, kin

        
    def add_mass(self, mass):
        self.mass += mass
    
####CORE

#Read input data
#file = 'Fonnbu_dhm.asc'
path = '/home/neuhauser/Documents/Milan/'
file = path + 'r_dem.tif'
release_file = path + 'r_src.tif'
#mass_out = 'Mass_exp1.asc'
vel_out =path + 'Velocity_v3.asc'
dem, header = io.read_raster(file)     
release, header = io.read_raster(release_file) 

#Initialize Output raster
mass_dist = np.zeros_like(dem)
vel_sqr_global = np.zeros_like(dem)
kin_energy = np.zeros_like(dem)

#for i in range(10):
#    dem[:,i] = dem[:,i] * (10-i)

start = time.time()

cell_list = []  
#startcell = Cell(500, 1, dem[500, 1], 0)
row_idx, col_idx = np.where(release == 1)
row_idx = row_idx[0]
col_idx = col_idx[0]
startcell = Cell(row_idx, col_idx, dem[row_idx, col_idx], 1, 1)
cell_list.append(startcell)
# ToDO: set ng_cells for array!!!
#ng_cells = startcell.rowindex-1:startcell.rowindex+2,startcell.colindex-1:startcell.colindex+2
# =============================================================================
# dem_ng = dem[startcell.rowindex-1:startcell.rowindex+2,startcell.colindex-1:startcell.colindex+2]
# velocity_sqr = startcell.calc_velocity(startcell, dem_ng)z i − 2z i+1
# for i in range(-1, 2):
#     for j in range(-1, 2):
#         if velocity_sqr[i+1, j+1] > vel_sqr[startcell.rowindex+i,startcell.colindex+j]:
#             vel_sqr[startcell.rowindex+i,startcell.colindex+j] = velocity_sqr[i+1, j+1]
# =============================================================================
        
            
#vel[ng_cells] = velocity

#mass_dist[startcell.rowindex, startcell.colindex] = 1
exist = False


for cell in cell_list:
    dem_ng = dem[cell.rowindex-1:cell.rowindex+2,cell.colindex-1:cell.colindex+2]
    velocity_sqr_local, kin = cell.calc_velocity(cell, dem_ng)
    row = []
    col = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            if velocity_sqr_local[i+1, j+1] > vel_sqr_global[cell.rowindex+i,cell.colindex+j]: # add velocity_sqr to global velocity_sqr raster
                vel_sqr_global[cell.rowindex+i,cell.colindex+j] = velocity_sqr_local[i+1, j+1]
                row.append(cell.rowindex+i)
                col.append(cell.colindex+j)

    for i in range(len(cell_list)): #Taking out multiple cells, ToDo: add velocity
         j = 0
         while j < len(row):
             if row[j] == cell_list[i].rowindex and col[j] == cell_list[i].colindex:
                 cell_list[i].add_mass(mass[j])
                 cell_list[i].velocity_sqr = max(cell_list[i].velocity_sqr, vel_sqr[row[j], col[j]])
                 if cell_list[i].velocity_sqr < vel_sqr_global[row[j], col[j]]:
                     cell_list[i].velocity_sqr = vel_sqr_global[row[j], col[j]]
                     
                 row = np.delete(row, j)
                 col = np.delete(col, j)
                 #mass = np.delete(mass, j)
             else:
                 j += 1

    for k in range(len(row)):
             cell_list.append(Cell(row[k], col[k], dem[row[k], col[k]], vel_sqr_global[row[k], col[k]]))
#
# """
#                 k = 0
#
#                 while k < len(cell_list): #CHeck for multiple cells
#                     bug = cell.rowindex+i
#                     bug1 = cell_list[k].rowindex
#                     bug2 = cell.colindex+j
#                     bug3 = cell_list[k].colindex
#                     if (cell.rowindex+i == cell_list[k].rowindex and cell.colindex+j == cell_list[k].colindex):
#                             #if cell_list[k].velocity_sqr < vel_sqr[cell.rowindex+i,cell.colindex+j]:
#                         cell_list[k].velocity_sqr = max(cell_list[k].velocity_sqr,vel_sqr[cell.rowindex+i,cell.colindex+j])
#                         k = len(cell_list)
#                         exist = True
#                     else:
#                         #cell_list.append(Cell(cell.rowindex+i, cell.colindex+j, dem[cell.rowindex+i, cell.colindex+j], velocity_sqr[i+1, j+1]))
#                         k +=1
#                 if not exist:
#                     cell_list.append(Cell(cell.rowindex + i, cell.colindex + j, dem[cell.rowindex + i, cell.colindex + j], velocity_sqr[i + 1, j + 1]))
# """
#


# =============================================================================
# row, col, mass = startcell.calc_distribution(dem[startcell.rowindex-1:startcell.rowindex+2,startcell.colindex-1:startcell.colindex+2])
# for k in range(len(row)): 
#     cell_list.append(Cell(row[k], col[k], dem[row[k], col[k]], mass[k]))
# 
# for cells in cell_list:
#     if not cells == startcell:
#         cells.calc_velocity(startcell)
#         vel[cells.rowindex, cells.colindex] = cells.velocity
#         mass_dist[cells.rowindex, cells.colindex] = cells.mass
#         if cells.velocity > 0:
#             row, col, mass = cells.calc_distribution(dem[cells.rowindex-1:cells.rowindex+2,cells.colindex-1:cells.colindex+2])
#             
#             for i in range(len(cell_list)): #Taking out multiple cells, ToDo: add velocity
#                 j = 0
#                 while j < len(row):
#                     if row[j] == cell_list[i].rowindex and col[j] == cell_list[i].colindex:
#                         cell_list[i].add_mass(mass[j])
#                         row = np.delete(row, j)
#                         col = np.delete(col, j)
#                         mass = np.delete(mass, j)
#                     else:
#                         j += 1
#                         
#             for k in range(len(row)):             
#                     cell_list.append(Cell(row[k], col[k], dem[row[k], col[k]], mass[k]))
# =============================================================================

end = time.time()            
print(end - start)
#ToDO: Implement if NoData Value is hit, or boarder of DEM            
# =============================================================================
#     k = 0    
#     while k < len(i): # Controll if cell is on boarder of DEM, when delete it from the list
#         if i[k] == np.size(dem, 0) - 1 or j[k] == np.size(dem, 1) - 1:
#             i = np.delete(i, k)
#             j = np.delete(j, k)
#         k += 1
# =============================================================================

    
raster_trans = rasterio.open(file)
crs = rasterio.crs.CRS.from_epsg(3912)
# =============================================================================
# new_dataset = rasterio.open(mass_out, 'w', driver='GTiff', height = mass_dist.shape[0], width = mass_dist.shape[1], count=1,  dtype = mass_dist.dtype, crs='+proj=latlong', transform=raster_trans.transform)
# new_dataset.write(mass_dist, 1)
# new_dataset.close()    
# =============================================================================
vel_sqr_global = np.sqrt(vel_sqr_global)
new_dataset = rasterio.open(vel_out, 'w', driver='GTiff', height = mass_dist.shape[0], width = mass_dist.shape[1], count=1,  dtype = mass_dist.dtype, crs=crs, transform=raster_trans.transform)
new_dataset.write(vel_sqr_global, 1)
new_dataset.close()    

