#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:14:39 2019

@author: Michael Neuhauser
This is the flow class
"""

import numpy as np


class Cell:
    
    def __init__(self, process,  rowindex, colindex, dem_ng, cellsize, susceptibility, elh, parent, alpha, exp, startcell):
        '''This class handles the spreading over the DEM!
        Depending on the process different alpha angles are used for energy dissipation.'''
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = dem_ng[1, 1]
        self.dem_ng = dem_ng
        self.cellsize = cellsize
        self.tan_beta = np.zeros_like(self.dem_ng)
        self.dist = np.zeros_like(self.dem_ng)
        self.persistence = np.zeros_like(self.dem_ng)
        self.p_fd = np.zeros_like(self.dem_ng)
        self.susceptibility = susceptibility
        self.elh = elh
        self.alpha = float(alpha)
        self.exp = int(exp)

        if process == 'Avalanche':
            self.p_threshold = 3 * 10 ** -4
            self.max_elh = 270  # maximum velocity this process can reach
        if process == 'Rockfall':
            self.p_threshold = 3 * 10 ** -4
            self.max_elh = 50  # maximum velocity this process can reach
        if process == 'Soil Slides':
            self.p_threshold = 3 * 10 ** -4
            self.max_elh = 12  # maximum velocity this process can reach
        self.parent = []
        if parent:
            self.parent.append(parent)
        
        if startcell:  # check, if start cell exist (start cell is release point)
            self.is_start = True  # set is_start to True
        else:            
            self.startcell = startcell  # give startcell to cell
            self.is_start = False  # set is_start to False
        
        self.calc_elh()
        self.calc_persistence()
        self.calc_tanbeta()

    def add_os(self, susceptibility):
        self.susceptibility += susceptibility
        
    def add_parent(self, parent):
        self.parent.append(parent)
        self.calc_persistence()
    
    def calc_elh(self):

        delta_e_pot = (self.dem_ng - self.altitude) * (-1)
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 0, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        tan_alpha = np.tan(np.deg2rad(self.alpha))
        e_friction = ds * self.cellsize * tan_alpha
        self.elh_neighbour = self.elh + delta_e_pot - e_friction
        self.elh_neighbour[self.elh_neighbour < 0] = 0
        self.elh_neighbour[self.elh_neighbour > self.max_elh] = self.max_elh
           
    def calc_tanbeta(self):  

        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 1, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        distance = ds * self.cellsize
        
        beta = np.arctan(((self.dem_ng - self.altitude) * (-1)) / distance) + 90
        self.tan_beta = np.tan(beta/2)

        self.tan_beta[self.elh_neighbour <= 0] = 0
        self.tan_beta[self.persistence <= 0] = 0
        self.tan_beta[1, 1] = 0
        if np.sum(self.tan_beta > 0):
            self.p_fd = self.tan_beta ** self.exp / np.sum(self.tan_beta ** self.exp)

    def calc_persistence(self):
        self.persistence = np.zeros_like(self.dem_ng)
        if self.is_start:
            self.persistence += 1
        elif self.parent[0].is_start:
            self.persistence += 1
        else:
            for parent in self.parent:
                dx = (parent.colindex - self.colindex) * -1 + 1
                dy = (parent.rowindex - self.rowindex) * -1 + 1
                maxweight = 1 * parent.elh
                if dx == 0 and dy == 0:
                    self.persistence[0, 0] += maxweight
                    self.persistence[1, 0] += 0.707 * maxweight
                    self.persistence[0, 1] += 0.707 * maxweight
                elif dx == 1 and dy == 0:
                    self.persistence[0, 1] += maxweight
                    self.persistence[0, 0] += 0.707 * maxweight
                    self.persistence[0, 2] += 0.707 * maxweight
                elif dx == 2 and dy == 0:
                    self.persistence[0, 2] += maxweight
                    self.persistence[0, 1] += 0.707 * maxweight
                    self.persistence[1, 2] += 0.707 * maxweight
                elif dx == 0 and dy == 1:
                    self.persistence[1, 0] += maxweight
                    self.persistence[2, 0] += 0.707 * maxweight
                    self.persistence[0, 0] += 0.707 * maxweight
                elif dx == 2 and dy == 1:
                    self.persistence[1, 2] += maxweight
                    self.persistence[0, 2] += 0.707 * maxweight
                    self.persistence[2, 2] += 0.707 * maxweight
                elif dx == 0 and dy == 2:
                    self.persistence[2, 0] += maxweight
                    self.persistence[1, 0] += 0.707 * maxweight
                    self.persistence[2, 1] += 0.707 * maxweight
                elif dx == 1 and dy == 2:
                    self.persistence[2, 1] += maxweight
                    self.persistence[2, 0] += 0.707 * maxweight
                    self.persistence[2, 2] += 0.707 * maxweight
                elif dx == 2 and dy == 2:
                    self.persistence[2, 2] += maxweight
                    self.persistence[2, 1] += 0.707 * maxweight
                    self.persistence[1, 2] += 0.707 * maxweight

            np.rot90(self.persistence, 2)
                    
    def calc_distribution(self):
        threshold = self.p_threshold
        if np.sum(self.p_fd) > 0:
            self.dist = self.persistence * self.p_fd / np.sum(self.persistence * self.p_fd) * self.susceptibility
        # This lines handle if a distribution to a neighbour cell is lower then the threshold, so we don´t lose
        # susceptibility.
        # The susceptibility of this cells will then spread equally to all neighbour cells
        count = ((0 < self.dist) & (self.dist < threshold)).sum()
        mass_to_distribute = np.sum(self.dist[self.dist < threshold])
        '''Checking if susceptibility is distributed to a field that isn't taking in account, when then distribute it to
         the other fields'''
        if mass_to_distribute > 0 and count > 0:
            self.dist[self.dist > threshold] += mass_to_distribute / count
            self.dist[self.dist < threshold] = 0
        if np.sum(self.dist) < self.susceptibility and count > 0:
            self.dist[self.dist > threshold] += (self.susceptibility - np.sum(self.dist))/count

        row_local, col_local = np.where(self.dist > threshold)
        
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.elh_neighbour[row_local, col_local]
