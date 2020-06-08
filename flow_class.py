#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:14:39 2019

@author: Michael Neuhauser
This is the graviclass
"""

import numpy as np


class Cell:
    
    def __init__(self, process,  rowindex, colindex, dem_ng, cellsize, mass, elh, forest, parent, alpha, exp, startcell):
        '''This class handles the spreading over the DEM!
        Depending on the process different alpha angles are used for energy dissipation.'''
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = dem_ng[1, 1]
        self.dem_ng = dem_ng
        self.cellsize = cellsize
        self.tan_beta = np.zeros_like(self.dem_ng)
        self.dist = np.zeros_like(self.dem_ng)
        self.direction = np.zeros_like(self.dem_ng)
        self.p_fd = np.zeros_like(self.dem_ng)
        self.forest = forest
        self.mass = mass
        self.kin_e = elh
        self.alpha = float(alpha)
        self.exp = int(exp)

        if process == 'Avalanche':
            #self.alpha = 25
            self.alpha_forest = 10
            #self.exp = 8
            self.mass_threshold = 3 * 10 ** -4
            self.max_elh = 270  # maximum velocity this process can reach
        if process == 'Rockfall':
            #self.alpha = 32
            self.alpha_forest = 0
            #self.exp = 75
            self.mass_threshold = 3 * 10 ** -4
            self.max_elh = 50  # maximum velocity this process can reach
        if process == 'Soil Slides':
            #self.alpha = 22
            self.alpha_forest = 0
            #self.exp = 75
            self.mass_threshold = 3 * 10 ** -4
            self.max_elh = 12  # maximum velocity this process can reach
        self.parent = []
        if parent:
            self.parent.append(parent)
        
        if startcell:  # check, if start cell exist (start cell is release point)
            self.is_start = True  # set is_start to True
        else:            
            self.startcell = startcell  # give startcell to cell
            self.is_start = False  # set is_start to False
        
        self.calc_kinetic_energy()
        self.calc_direction()
        self.calc_tanbeta()

    def add_mass(self, mass):
        self.mass += mass
        
    def add_parent(self, parent):
        self.parent.append(parent)
        self.calc_direction()
    
    def calc_kinetic_energy(self):

        delta_e_kin_pot = (self.dem_ng - self.altitude) * (-1)
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 0, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        max_friction = self.alpha_forest * self.forest
        if self.kin_e < 45:
            alpha_calc = self.alpha + max(0, - self.kin_e * (max_friction / 45) + max_friction)
        else:
            alpha_calc = 25
        # print(alpha_calc)
        tan_alpha = np.tan(np.deg2rad(alpha_calc))  # increased friction due to forest scaled with forest value (forest)
        e_friction = ds * self.cellsize * tan_alpha
        self.kin_energy_neighbour = self.kin_e + delta_e_kin_pot - e_friction
        self.kin_energy_neighbour[self.kin_energy_neighbour < 0] = 0
        self.kin_energy_neighbour[self.kin_energy_neighbour > self.max_elh] = self.max_elh
                    
           
    def calc_tanbeta(self):  

        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 1, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        distance = ds * self.cellsize
        
        beta = np.arctan(((self.dem_ng - self.altitude) * (-1)) / distance) + 90
        self.tan_beta = np.tan(beta/2)

        self.tan_beta[self.kin_energy_neighbour <= 0] = 0
        self.tan_beta[self.direction <= 0] = 0
        self.tan_beta[1, 1] = 0
        if np.sum(self.tan_beta > 0):
            self.p_fd = self.tan_beta ** self.exp / np.sum(self.tan_beta ** self.exp)

    def calc_direction(self):
        self.direction = np.zeros_like(self.dem_ng)
        if self.is_start:
            self.direction += 1
        elif self.parent[0].is_start:
            self.direction += 1
        else:
            for parent in self.parent:
                dx = (parent.colindex - self.colindex) * -1 + 1
                dy = (parent.rowindex - self.rowindex) * -1 + 1
                maxweight = 1 * parent.kin_e
                if dx == 0 and dy == 0:
                    self.direction[0, 0] += maxweight
                    self.direction[1, 0] += 0.707 * maxweight
                    self.direction[0, 1] += 0.707 * maxweight
                elif dx == 1 and dy == 0:
                    self.direction[0, 1] += maxweight
                    self.direction[0, 0] += 0.707 * maxweight
                    self.direction[0, 2] += 0.707 * maxweight
                elif dx == 2 and dy == 0:
                    self.direction[0, 2] += maxweight
                    self.direction[0, 1] += 0.707 * maxweight
                    self.direction[1, 2] += 0.707 * maxweight
                elif dx == 0 and dy == 1:
                    self.direction[1, 0] += maxweight
                    self.direction[2, 0] += 0.707 * maxweight
                    self.direction[0, 0] += 0.707 * maxweight
                elif dx == 2 and dy == 1:
                    self.direction[1, 2] += maxweight
                    self.direction[0, 2] += 0.707 * maxweight
                    self.direction[2, 2] += 0.707 * maxweight
                elif dx == 0 and dy == 2:
                    self.direction[2, 0] += maxweight
                    self.direction[1, 0] += 0.707 * maxweight
                    self.direction[2, 1] += 0.707 * maxweight
                elif dx == 1 and dy == 2:
                    self.direction[2, 1] += maxweight
                    self.direction[2, 0] += 0.707 * maxweight
                    self.direction[2, 2] += 0.707 * maxweight
                elif dx == 2 and dy == 2:
                    self.direction[2, 2] += maxweight
                    self.direction[2, 1] += 0.707 * maxweight
                    self.direction[1, 2] += 0.707 * maxweight

            np.rot90(self.direction,2)
                    
    def calc_distribution(self):
        threshold = self.mass_threshold
        if np.sum(self.p_fd > 0):
            self.dist = self.direction * self.p_fd / np.sum(self.direction * self.p_fd) * self.mass
        # This lines handle if a distribution to a neighbour cell is lower then the threshold, so we don´t lose "mass"
        # The mass of this cells will then be spreaded equally to all neighbour cells
        count = ((0 < self.dist) & (self.dist < threshold)).sum()
        mass_to_distribute = np.sum(self.dist[self.dist < threshold])
        '''Checking if mass is distributed to a field that isnt taking in account, when then distribute it to the other 
        fields'''
        if mass_to_distribute > 0 and count > 0:
            self.dist[self.dist > threshold] += mass_to_distribute / count
            self.dist[self.dist < threshold] = 0
        if np.sum(self.dist) < self.mass and count > 0:
            self.dist[self.dist > threshold] += (self.mass - np.sum(self.dist))/count

        row_local, col_local = np.where(self.dist > threshold)
        
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.kin_energy_neighbour[row_local, col_local]
