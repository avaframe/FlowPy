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

class cell():
    
    def __init__(self, rowindex, colindex, altitude, mass):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = altitude
        self.mass = mass
        self.velocity = 0
        self.cellsize = 10
        self.exp = 1
        self.alpha = 50
                       
    def calc_velocity(self, startcell):
        dx = np.abs(startcell.colindex - self.colindex) * self.cellsize
        dy = np.abs(startcell.rowindex - self.rowindex) * self.cellsize
        dh = np.sqrt(dx**2 + dy**2)
        temp = startcell.altitude - self.altitude - dh * np.tan(np.deg2rad(self.alpha))
        if temp > 0:
            self.velocity = np.sqrt(2*9.81*temp)
        else:
            self.velocity = 0
            
        return self.velocity
        
    def add_mass(self, mass):
        self.mass += mass
        
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
                
    def calc_distribution(self, dem_ng):      
        self.calc_tanbeta(dem_ng)
        self.mass_dist = np.zeros((3,3))
        exp = self.exp
        mass_threshold = 0.0001

        self.tan_beta[self.tan_beta > 0] = 0
        for i in range(3):
            for j in range(3):
                self.mass_dist[i, j] = self.mass * self.tan_beta[i, j] ** exp / np.sum(self.tan_beta ** exp)
        out_row_local, out_col_local = np.where(self.mass_dist > mass_threshold)
        out_row = self.rowindex - 1 + out_row_local
        out_col = self.colindex - 1 + out_col_local
        mass_list = self.mass_dist[out_row_local, out_col_local]
        
        return out_row, out_col, mass_list
    
    
def read_header(input_file):

    raster = rasterio.open(input_file)
    if raster is None:
        print('Unable to open {}'.format(input_file))
        sys.exit(1)

    header = {}
    header['ncols'] = raster.width
    header['nrows'] = raster.height
    header['xllcorner'] = (raster.transform * (0, 0))[0]
    header['yllcorner'] = (raster.transform * (0, raster.height))[1]
    header['cellsize'] = raster.transform[0]
    header['noDataValue'] = raster.nodata
    return header


def read_raster(input_file):

    header = read_header(input_file)
    raster = rasterio.open(input_file)
    my_array = raster.read(1)

    return my_array, header


file = 'Fonnbu_dhm.asc'
mass_out = 'Mass_exp8.asc'
vel_out = 'Velocity.asc'
dem, header = read_raster(file)   
#dem = np.ones((100, 100), dtype=np.int)
mass_dist = np.zeros_like(dem)
vel = np.zeros_like(dem)

#for i in range(10):
#    dem[:,i] = dem[:,i] * (10-i)

start = time.time()

cell_list = []  
startcell = cell(500, 1, dem[500, 1], 1)
#mass_dist[startcell.rowindex, startcell.colindex] = 1

cell_list.append(startcell)
row, col, mass = startcell.calc_distribution(dem[startcell.rowindex-1:startcell.rowindex+2,startcell.colindex-1:startcell.colindex+2])
for k in range(len(row)): 
    cell_list.append(cell(row[k], col[k], dem[row[k], col[k]], mass[k]))

for cells in cell_list:
    if not cells == startcell:
        vel[cells.rowindex, cells.colindex] = cells.calc_velocity(startcell)
        if cells.velocity > 0:
            print(cells.velocity)
            row, col, mass = cells.calc_distribution(dem[cells.rowindex-1:cells.rowindex+2,cells.colindex-1:cells.colindex+2])
            
            for i in range(len(cell_list)): #Taking out multiple cells, ToDo: add velocity
                j = 0
                while j < len(row):
                    if row[j] == cell_list[i].rowindex and col[j] == cell_list[i].colindex:
                        cell_list[i].add_mass(mass[j])
                        row = np.delete(row, j)
                        col = np.delete(col, j)
                        mass = np.delete(mass, j)
                    else:
                        j += 1
                        
            for k in range(len(row)):             
                    cell_list.append(cell(row[k], col[k], dem[row[k], col[k]], mass[k]))
                
            mass_dist[cells.rowindex, cells.colindex] = cells.mass
        else:
            break
    #vel[cells.rowindex, cells.colindex] = cells.calc_velocity(startcell)

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
new_dataset = rasterio.open(mass_out, 'w', driver='GTiff', height = mass_dist.shape[0], width = mass_dist.shape[1], count=1,  dtype = mass_dist.dtype, crs='+proj=latlong', transform=raster_trans.transform)
new_dataset.write(mass_dist, 1)
new_dataset.close()    

new_dataset = rasterio.open(vel_out, 'w', driver='GTiff', height = mass_dist.shape[0], width = mass_dist.shape[1], count=1,  dtype = mass_dist.dtype, crs='+proj=latlong', transform=raster_trans.transform)
new_dataset.write(vel, 1)
new_dataset.close()    

