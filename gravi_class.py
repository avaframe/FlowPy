#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:14:39 2019

@author: neuhauser
This is the graviclass
"""

import numpy as np


class Cell():
    
    def __init__(self, rowindex, colindex, dem_ng, cellsize, mass, kin_e, forest, parent, startcell):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = dem_ng[1, 1]
        self.dem_ng = dem_ng
        self.cellsize = cellsize
        self.tan_beta = np.zeros_like(self.dem_ng)
        self.dist = np.zeros_like(self.dem_ng)
        self.direction = np.zeros_like(self.dem_ng)
        self.p_fd = np.zeros_like(self.dem_ng)
        self.alpha = 25
        self.alpha_forest = 5
        self.exp = 8 
        self.mass_threshold = 3*10**-4
        self.forest = forest
        self.mass = mass
        self.kin_e = kin_e
        self.parent = []
        if parent:
            self.parent.append(parent)
        
        if startcell == True: #check, if start cell exist (start cell is release point)
            self.is_start = True # set is_satrt to True
        else:            
            self.startcell = startcell # set is_satrt to True
            self.is_start = False # set is_satrt to False        
        
        self.calc_kinetic_energy()
        self.calc_direction()
        #self.calc_global_direction()
        self.calc_tanbeta()
        
       
    def calc_kinetic_energy(self):
        delta_e_kin_pot = (self.dem_ng - self.altitude) * (-1)
        ds = np.array([[np.sqrt(2),1,np.sqrt(2)],[1,0,1],[np.sqrt(2),1,np.sqrt(2)]])
        tan_alpha = np.tan(np.deg2rad(self.alpha + self.forest * self.alpha_forest))  # increased friction due to forest scaled with forest value (forest)
        e_friction = ds * self.cellsize * tan_alpha
        self.kin_energy_neighbour = self.kin_e + delta_e_kin_pot - e_friction
        self.kin_energy_neighbour[self.kin_energy_neighbour < 0] = 0
                    
    def add_mass(self, mass):
        self.mass += mass
        
    def add_parent(self, parent):
        self.parent.append(parent)
        self.calc_direction()
           
    def calc_tanbeta(self):  
        exp = self.exp
        #snowdepth = 1
        #density = 100
        #dh = self.kin_e/(9.81*self.mass*self.cellsize**2 * snowdepth * density) # Calculate the remaining Energyheight with a kind of mass.... 
# =============================================================================
#         if self.is_start:
#             dh = 0
#         else:
#             dx = (self.startcell.colindex - self.colindex) * self.cellsize
#             dy = (self.startcell.rowindex - self.rowindex) * self.cellsize
#             ds = np.sqrt(dx**2 + dy**2)
#             #dh = (self.startcell.altitude - self.altitude - ds * np.tan(np.deg2rad(self.alpha)))
#             #dh = self.kin_e / 9.81
#             dh = 10
# =============================================================================

        ds = np.array([[np.sqrt(2),1,np.sqrt(2)],[1,0,1],[np.sqrt(2),1,np.sqrt(2)]])
        distance = ds * self.cellsize
        dh = 1
        self.tan_beta = ((self.dem_ng - (self.altitude + dh)) * (-1)) / distance
        #self.tan_beta = np.tan((beta+90)/2)

        #self.tan_beta[self.tan_beta < 0] = 0
        #self.tan_beta = abs(self.tan_beta)
        self.tan_beta[self.kin_energy_neighbour <= 0] = 0
        self.tan_beta[self.direction <= 0] = 0
        self.tan_beta[1,1] = 0
        if np.sum(self.tan_beta > 0):
            self.p_fd = self.tan_beta ** exp/ np.sum(self.tan_beta ** exp)
        
    def calc_direction(self):
        if self.is_start:
            self.direction += 1
        #elif self.parent[0].is_start:
            #self.direction += 1
        else:
            for parent in self.parent:
                dx = (parent.colindex - self.colindex) * -1 +1
                dy = (parent.rowindex - self.rowindex) * -1 +1
                maxweight = 1
                if dx == 0 and dy == 0:
                    self.direction[0, 0] += maxweight
                    self.direction[1, 0] += 0.707
                    self.direction[0, 1] += 0.707
                elif dx == 1 and dy == 0:
                    self.direction[0, 1] += maxweight
                    self.direction[0, 0] += 0.707
                    self.direction[0, 2] += 0.707
                elif dx == 2 and dy == 0:
                    self.direction[0, 2] += maxweight
                    self.direction[0, 1] += 0.707
                    self.direction[1, 2] += 0.707
                elif dx == 0 and dy == 1:
                    self.direction[1, 0] += maxweight
                    self.direction[2, 0] += 0.707
                    self.direction[0, 0] += 0.707
                elif dx == 2 and dy == 1:
                    self.direction[1, 2] += maxweight
                    self.direction[0, 2] += 0.707
                    self.direction[2, 2] += 0.707
                elif dx == 0 and dy == 2:
                    self.direction[2, 0] += maxweight
                    self.direction[1, 0] += 0.707
                    self.direction[2, 1] += 0.707
                elif dx == 1 and dy == 2:
                    self.direction[2, 1] += maxweight
                    self.direction[2, 0] += 0.707
                    self.direction[2, 2] += 0.707
                elif dx == 2 and dy == 2:
                    self.direction[2, 2] += maxweight
                    self.direction[2, 1] += 0.707
                    self.direction[1, 2] += 0.707

            np.rot90(self.direction,2)
            
    def calc_global_direction(self):
        
        if self.is_start:
            self.global_dir = np.ones((3,3))
        else:
            dx = (self.colindex - self.startcell.colindex) * self.cellsize # x component of avalanche flow direction, global
            dy = (self.rowindex - self.startcell.rowindex) * self.cellsize # y component of avalanche flow direction, global
            avi_direction = np.rad2deg(np.arctan2(dy, dx)) * -1 # avalanche direction in degrees, global
            #neighbour_direction = np.linspace(0.0, 315, 8)  # direction in deg to all cells
            
            # Setting the direction for the Center Cell in the opositre direction for the avalanche, so it´s not calculated
            neighbour_direction = np.array([[225, 270, 315], [180, np.nan, 0], [225, 270 , 315]])  # direction in deg to all cells, local
            delta_deg = neighbour_direction - avi_direction  # difference between ava direction and the neighbor cells
            self.global_dir = np.cos(np.deg2rad(delta_deg))
            #self.global_dir *= 0.3# direction_projection[1] = NE, direction_projection[2] = N ..., global influence
            self.global_dir[self.global_dir < 0.1] = 0
            self.global_dir[1, 1] = 0
           
                    
    def calc_distribution(self):
        threshold = self.mass_threshold
        if np.sum(self.p_fd > 0):
            self.dist = self.direction * self.p_fd / np.sum(self.direction * self.p_fd) * self.mass
            #self.dist = self.global_dir * self.p_fd / np.sum(self.global_dir * self.p_fd) * self.mass
            #self.dist = (self.direction*self.kin_e + self.tan_beta)/np.sum(self.direction*self.kin_e + self.tan_beta)*self.mass
        count = ((0 < self.dist) & (self.dist < threshold)).sum()
        mass_to_distribute = np.sum(self.dist[self.dist < threshold])
        if mass_to_distribute > 0 and count > 0:
            self.dist[self.dist > threshold] += mass_to_distribute / count
            self.dist[self.dist < threshold] = 0
        if np.sum(self.dist) < self.mass and count > 0:
            self.dist[self.dist > threshold] += (self.mass - np.sum(self.dist))/count
            #print('Mass Loss' , np.sum(self.dist) - self.mass)
        row_local, col_local = np.where(self.dist > threshold)  # Zellen die nicht im threshold liegen müssen ihre masse auf die anderen verteilen!
        
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.kin_energy_neighbour[row_local, col_local]
 