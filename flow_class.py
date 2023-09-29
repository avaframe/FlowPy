#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:14:39 2019

This is the flow class

    Copyright (C) <2020>  <Michael Neuhauser>
    Michael.Neuhauser@bfw.gv.at

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
import math


class Cell:
    #def __init__(self, rowindex, colindex, dem_ng, forest, cellsize, flux, z_delta, parent, alpha, exp, flux_threshold, max_z_delta, startcell):
    #PAULA
    def __init__(self, rowindex, colindex, dem_ng, forest, cellsize, flux, z_delta, parent, alpha, exp, flux_threshold, max_z_delta, startcell, co_e, co_m):
    #ende paula
        '''This class handles the spreading over the DEM!
        Depending on the process different alpha angles are used for energy dissipation.'''
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = dem_ng[1, 1]
        self.dem_ng = dem_ng
        self.forest = forest
        self.cellsize = cellsize
        self.tan_beta = np.zeros_like(self.dem_ng)
        self.dist = np.zeros_like(self.dem_ng)
        self.persistence = np.zeros_like(self.dem_ng)
        self.r_t = np.zeros_like(self.dem_ng)
        self.no_flow = np.ones_like(self.dem_ng)
        self.flux = flux
        self.z_delta = z_delta
        #JT
        self.altitutde_diff = 0
        #ende
        #PAULA
        self.flow_energy = 0
        self.dist_energy = np.zeros_like(self.dem_ng)
        self.coe_energy_avg = 0 # average value of flow energy, if cell is in center of line
        self.com_energy_avg = 0
        self.coe_flux_avg = 0 # average value of flux, if cell is in center of line
        self.com_flux_avg = 0
        self.co_e = co_e
        self.co_m = co_m #size: dem, is 1 if cell is in center of
        self.horizontal_diff = 0
        self.coe_ds = 0 #size: dem, distance to next center of cell
        self.row_idx_coe = 0 #row index of next center of cell
        self.col_idx_coe = 0 #col index of next center of cell
        self.row_idx_com = 0 #row index of next center of cell
        self.col_idx_com = 0 #col index of next center of cell
        self.s_com = 0 # s calculated weighted via center of mass
        self.s_coe = 0
        self.iteration = 0
        #ende paula
        self.alpha = float(alpha)
        self.exp = int(exp)
        self.max_z_delta = float(max_z_delta)
        self.flux_threshold = float(flux_threshold)
        self.min_distance = 0
        self.max_distance = 0
        self.min_gamma = 0
        self.max_gamma = 0
        self.sl_gamma = 0    
        # Parameters for Forest Friction, right now use it just for avalanches   
        #self.alpha_forest = 10  # Max added friction angel
        self.max_added_friction_forest = 10 # degrees added to friction angle
        self.min_added_friction_forest = 2 # minimium effect forested terrain can have
        self.no_friction_effect_v = 30 # velocity shared for friction and detrainment methods
        self.max_added_detrainment_forest = 0.0003 #
        self.min_added_detrainment_forest = 0#0.00001
        self.no_detrainmnet_effect_v = 30 #

        if type(startcell) == bool:  # check, if start cell exist (start cell is release point)
            self.is_start = True  # set is_start to True
        else:            
            self.startcell = startcell  # give startcell to cell
            self.is_start = False  # s et is_start to False

        self.parent = []
        if type(parent) == Cell:
            self.parent.append(parent)
            
    def add_os(self, flux):
        self.flux += flux

    def forest_detrainment(self ):
        """
        linear decrease of forest effect with regard to alpha increase and kinetic energy height
        This is the detrainment routine for forest. It should reduce the routing flux of the avalanche.
        self.max_added_detrainment_forest = .001
        self.min_added_detrainment_forest = 0
        self.no_detrainmnet_effect_v = 45
        """
        no_detrainmnet_effect_zdelta = self.no_detrainmnet_effect_v**(2)/(np.sqrt(2) * 9.8)  # change veloctiy into kinetic energy line hight (z_delta) 9.8 = gravity. derrived from 1/2mv^2 = mgh
        rest = self.max_added_detrainment_forest * self.forest # detrainment effect scalled to forest, should be zero for non-forested area
        slope = (rest - self.min_added_detrainment_forest)/(0-no_detrainmnet_effect_zdelta) # rise over run (should be negative slope)
        self.detrainment  = max(self.min_added_detrainment_forest, slope * self.z_delta + rest) # y = mx + b, shere z_delta is the x


    def add_parent(self, parent):
        self.parent.append(parent)

    def calc_altitude_diff(self):
        #JT
        self.altitutde_diff = self.startcell.altitude - self.altitude
        #ende

    def calc_fp_travelangle(self):
        dist_min = []
        dh = self.startcell.altitude - self.altitude
        self.altitutde_diff = dh
        for parent in self.parent:
            dx = abs(parent.colindex - self.colindex)
            dy = abs(parent.rowindex - self.rowindex)
            dist_min.append(math.sqrt(dx ** 2 + dy ** 2) * self.cellsize + parent.min_distance)
        self.min_distance = np.amin(dist_min)
        self.max_gamma = np.rad2deg(np.arctan(dh / self.min_distance))

    def calc_sl_travelangle(self):
        dx = abs(self.startcell.colindex - self.colindex)
        dy = abs(self.startcell.rowindex - self.rowindex)
        dh = self.startcell.altitude - self.altitude

        ds = math.sqrt(dx ** 2 + dy ** 2) * self.cellsize
        #PAULA
        self.horizontal_diff = ds
        #ende paula
        self.sl_gamma = np.rad2deg(np.arctan(dh / ds))

        
    def calc_theta(self):
        #new (Paula): calculate slope
        dz_dy, dz_dx = np.gradient(self.dem_ng, self.cellsize, self.cellsize)
        slope_rad = np.arctan(np.sqrt(dz_dx ** 2 + dz_dy ** 2))
        return slope_rad

    def calc_Voellmy_friction(self):
        #new (Paula): calculate turbulence term
        #assume constants
        g = 9.81 # m s-2
        rho = 200 # kg m-3
        
        muVoellmy = 0.155
        xsiVoellmy = 1000. #400-4000 (https://doi.org/10.3189/2015JoG14J168)  # 4000. (avaframe) #m/s²
        
        u = np.sqrt(self.z_delta * 2 * g)
        # mass = self.flux * cellsize**2 * h0 * rho
        # h = mass / cellsize**2 / rho
        h = self.flux * 10000 / self.cellsize**2
        V = self.cellsize ** 2 * h # m³ Volume
        
        theta = self.calc_theta()
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 0, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        
        tan_alpha_turb_voellmy = u * u * ds * self.cellsize / np.cos(theta) / h / xsiVoellmy
        return tan_alpha_turb_voellmy

    def calc_z_delta(self):
        self.z_delta_neighbour = np.zeros((3, 3))
        self.z_gamma = self.altitude - self.dem_ng
        no_friction_effect_zdelta = self.no_friction_effect_v**(2) / (np.sqrt(2) * 9.8) # change veloctiy into energy line hight (z_delta) 9.8 = gravity. derrived from 1/2mv^2 = mgh
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 0, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        ## Calculation for Forest Friction leads to new alpha_calc
        if self.forest > 0:
            if self.z_delta < no_friction_effect_zdelta: # no min_added forest values becuase of this line
                rest = self.max_added_friction_forest * self.forest  # friction at rest v=0 would be applied to start cells
                slope = (rest - self.min_added_friction_forest) / (0 - no_friction_effect_zdelta)  # rise over run
                friction = max(self.min_added_friction_forest,
                                   slope * self.z_delta + rest)  # y = mx + b, shere z_delta is the x
                alpha_calc = self.alpha + max(0, friction)
            else:
                alpha_calc = self.alpha + self.min_added_friction_forest
        else:
            alpha_calc = self.alpha
        # Normal calculation
        tan_alpha = np.tan(np.deg2rad(alpha_calc))
        self.z_alpha = ds * self.cellsize * tan_alpha
        
        #new (Paula): calculate friction including turbulence term (Voellmy)
        muVoellmy = 0.155
        #self.z_alpha = muVoellmy * ds * self.cellsize + self.calc_Voellmy_friction()
        #self.z_alpha += self.calc_Voellmy_friction()
        # END Paula
        
        
        self.z_delta_neighbour = self.z_delta + self.z_gamma - self.z_alpha
        self.z_delta_neighbour[self.z_delta_neighbour < 0] = 0
        self.z_delta_neighbour[self.z_delta_neighbour > self.max_z_delta] = self.max_z_delta
           
    def calc_tanbeta(self):
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 1, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        distance = ds * self.cellsize
        
        beta = np.arctan((self.altitude - self.dem_ng) / distance) + np.deg2rad(90)
        self.tan_beta = np.tan(beta/2)

        self.tan_beta[self.z_delta_neighbour <= 0] = 0
        self.tan_beta[self.persistence <= 0] = 0
        self.tan_beta[1, 1] = 0
        if abs(np.sum(self.tan_beta)) > 0:
            self.r_t = self.tan_beta ** self.exp / np.sum(self.tan_beta ** self.exp)
            
    def calc_flow_energy(self):
        ##NEW PAULA
    	self.flow_energy = self.flux * self.z_delta / 2

    def calc_iteration_step(self):
        #NEW PAULA
        if self.is_start == False:
            for parent in self.parent:
                self.iteration = parent.iteration + 1

    def calc_centerof(self, parameter):
        #NEW PAULA
        y = [[-1,-1,-1],[0,0,0],[1,1,1]]
        x = [[-1,0,1],[-1,0,1],[-1,0,1]]
        s = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 0, 1], [np.sqrt(2), 1, np.sqrt(2)]]) * self.cellsize
        if parameter == 'e':
            distr = self.dist_energy
            cell_value = self.flow_energy
        elif parameter == 'm':
            distr = self.dist
            cell_value = self.flux

        if cell_value > 0:
            tot = cell_value
        else: 
            tot = sum(sum(distr))
        x_co = 1/tot * sum(sum(x*distr))
        y_co = 1/tot * sum(sum(y*distr))
        row_idx_co = self.rowindex + round(y_co)
        col_idx_co = self.colindex + round(x_co)
        setattr(self,f's_co{parameter}', 1/tot * sum(sum(s*distr)))
        setattr(self, f'row_idx_co{parameter}', row_idx_co)
        setattr(self, f'col_idx_co{parameter}', col_idx_co)
        # self.s_co = 1/tot * sum(sum(s*distr))
        # self.row_idx_co = self.rowindex + round(y_co)
        # self.col_idx_co = self.colindex + round(x_co)
        try:
            setattr(self,f'co{parameter}_energy_avg',1/np.count_nonzero(self.dist_energy) * sum(sum(self.dist_energy)))
            #self.co_energy_avg = 1/np.count_nonzero(distr) * sum(sum(self.dist_energy))
            setattr(self,f'co{parameter}_flux_avg',1/np.count_nonzero(self.dist) * sum(sum(self.dist)))
        except: 
            setattr(self,f'co{parameter}_energy_avg',0)
            setattr(self,f'co{parameter}_flux_avg',0)
        if parameter == 'e':
            self.co_e[row_idx_co,col_idx_co] = 1
        elif parameter == 'm':
            self.co_m[row_idx_co,col_idx_co] = 1
        
        # Distance to next cell along center of ...
        dx = abs(self.colindex - col_idx_co)
        dy = abs(self.rowindex - row_idx_co)
        #self.coe_ds = math.sqrt(dx ** 2 + dy ** 2) * self.cellsize
        setattr(self, f'co{parameter}_ds', math.sqrt(dx ** 2 + dy ** 2) * self.cellsize)
        

    def calc_persistence(self):
        self.persistence = np.zeros_like(self.dem_ng)
        if self.is_start:
            self.persistence += 1
        elif self.parent[0].is_start:
            self.persistence += 1
        else:
            for parent in self.parent:
                dx = (parent.colindex - self.colindex) 
                dy = (parent.rowindex - self.rowindex)

                self.no_flow[dy + 1,dx + 1] = 0  # 3x3 Matrix of ones, every parent gets a 0, so no flow to a parent field.
                
                maxweight = parent.z_delta
                # Old Calculation
                if dx == -1:
                    if dy == -1:
                        self.persistence[2, 2] += maxweight
                        self.persistence[2, 1] += 0.707 * maxweight
                        self.persistence[1, 2] += 0.707 * maxweight
                    if dy == 0:
                        self.persistence[1, 2] += maxweight
                        self.persistence[2, 2] += 0.707 * maxweight
                        self.persistence[0, 2] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[0, 2] += maxweight
                        self.persistence[0, 1] += 0.707 * maxweight
                        self.persistence[1, 2] += 0.707 * maxweight

                if dx == 0:
                    if dy == -1:
                        self.persistence[2, 1] += maxweight
                        self.persistence[2, 0] += 0.707 * maxweight
                        self.persistence[2, 2] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[0, 1] += maxweight
                        self.persistence[0, 0] += 0.707 * maxweight
                        self.persistence[0, 2] += 0.707 * maxweight

                if dx == 1:
                    if dy == -1:
                        self.persistence[2, 0] += maxweight
                        self.persistence[1, 0] += 0.707 * maxweight
                        self.persistence[2, 1] += 0.707 * maxweight
                    if dy == 0:
                        self.persistence[1, 0] += maxweight
                        self.persistence[0, 0] += 0.707 * maxweight
                        self.persistence[2, 0] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[0, 0] += maxweight
                        self.persistence[0, 1] += 0.707 * maxweight
                        self.persistence[1, 0] += 0.707 * maxweight
                        
# =============================================================================
#                 # New Calculation:
#                 theta_child = np.array([[np.pi*5/4, np.pi*3/2 , np.pi*7/4], [np.pi, 0, 0], [np.pi*3/4, np.pi/2 , np.pi/4]])
#                 theta_parent = (np.arctan2(dy, dx))
#                 
#                 pers1 = theta_parent - theta_child - np.pi
#                 pers = np.zeros((3,3))
#                 
#                 for idx, element in np.ndenumerate(pers1):
#                     pers[idx] = max(0, np.cos(element))
#                     if pers[idx] < 2*np.finfo(np.float64).eps:
#                         pers[idx] = 0
#                 pers[1, 1] = 0
#                 self.persistence += pers * maxweight
# =============================================================================

    def calc_distribution(self):

        self.calc_iteration_step()
        self.calc_z_delta()
        self.calc_persistence()
        self.persistence *= self.no_flow
        self.calc_tanbeta()
        #print(self.persistence)
        self.forest_detrainment()
        

        if not self.is_start:
            self.calc_fp_travelangle()
            self.calc_sl_travelangle()
            #JT
            self.calc_altitude_diff()
            #ende

        self.flux = max(0.0003, self.flux - self.detrainment) # here we subtract the detrainment from the flux before moving flux to new cells.
        
        #PAULA
        self.calc_flow_energy()
     
        #ende Paula

        threshold = self.flux_threshold
        if np.sum(self.r_t) > 0: # if there is routing flux
            self.dist = (self.persistence * self.r_t) / np.sum(self.persistence * self.r_t) * self.flux
            #detrainment =  forest_scale(self, max_detrainment, min_detrainment, FSI, max_FE_velocity, z_delta )
            #self.dist = self.dist - self.detrainment # todo max needed to keep it always positive
            #todo make sure that I am only removing detrained mass once not every spatial iteration that includes a cell as a neighbor
            #print(self.detrainment, "detrainment" , self.dist, "mass")
        # This lines handle if a distribution to a neighbour cell is lower then the threshold, so we don´t lose
        # flux.
        # The flux of this cells will then spread equally to all neighbour cells
        count = ((0 < self.dist) & (self.dist < threshold)).sum()
        mass_to_distribute = np.sum(self.dist[self.dist < threshold])
        '''Checking if flux is distributed to a field that isn't taking in account, when then distribute it equally to
         the other fields'''
        if mass_to_distribute > 0 and count > 0:
            self.dist[self.dist > threshold] += mass_to_distribute / count
            self.dist[self.dist < threshold] = 0
        if np.sum(self.dist) < self.flux and count > 0:
            self.dist[self.dist > threshold] += (self.flux - np.sum(self.dist))/count
            
        self.dist_energy = self.dist * self.z_delta_neighbour / 2

        row_local, col_local = np.where(self.dist > threshold)
        #PAULA
        if self.co_e[self.rowindex, self.colindex] == 1:
            self.calc_centerof(parameter = 'e')
            #print(self.colindex, self.rowindex, self.dist_energy)
        if self.co_m[self.rowindex, self.colindex] == 1:
            self.calc_centerof(parameter = 'm')
            #print(self.colindex, self.rowindex, self.dist)
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.z_delta_neighbour[row_local, col_local], self.co_e, self.co_m
        #return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.z_delta_neighbour[row_local, col_local], self.rowindex - 1 + row_max, self.colindex - 1 + col_max
         #ende Paula
        #return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.z_delta_neighbour[row_local, col_local]
