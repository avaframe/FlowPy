#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
last modified 2023-07-14 A.H. (andreas.huber@bfw.gv.at)

This is an attempt to implement the original flow_clas.py Cell class
as a Cython Extension Module



########################################################################

Original code by: 
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
########################################################################    

"""

import cython
import numpy as np
cimport numpy as np
from libc cimport math as cMath

'''
some constants are pre-computed here, in order to avoid
unnecessary function calls in the code
'''
cdef double SQRT2 = cMath.sqrt(2.)
cdef double PI = cMath.M_PI
cdef double HALFPI = cMath.M_PI_2
cdef double QUATERPI = cMath.M_PI_4

cdef class Cell:
    """
    first try for a cython extension-type implementation of flow_class.py (based on master-branch 2023-07-10)
    This class handles spreading of the process on the level of the 3x3 neighbourhood, i.e. the 'cell level'    
    """
    
    #----Variables related to DEM and positioning within DEM--------#
    cdef public long rowindex, colindex #index of row and columnt of the cell in the raster-domain
    cdef public float altitude #elevation of central cell - from DEM
    cdef public float cellsize #cellsize of the computational domain/input DEM
    
    #----Flow-Py Model Parameters +- defined by user input--------#
    cdef public float alpha   #local alpha angle on the cell - either globally defined or locally adjusted (e.g. forest vs. non forest)
    cdef public float exp     #local spreading coefficient - global model parameter between 1 and inf
    cdef public float flux_threshold #minimum flux required for propagation - gobal model parameter
    cdef public float max_z_delta #limitation of kinetic energy / energy line height - proxy for max. velocity/velocity limitation
    
    #----Flow-Py Model Parameters +- Forest extension Chris D'Amboise--------#
    cdef public float FSI #'Forest Structure Index' between 0 and 1 used for scaling forest effects 0 -no forest 1 - ultra dense forest
    cdef public float min_afForest, max_afForest #min and max additional friction (in degrees) in forested areas
    cdef public float min_adForest, max_adForest #min and max detrainment values/pseudo Volumes in forested areas
    cdef public float maxEffectVForest_F #above this value (m/s) Forest friction effect on flow behavior assumed to be negligible
    cdef public float maxEffectVForest_D #above this value (m/s) Forest Detraiment effect on flow behavior assumed to be negligible 
    
    #-----helper variables/result variables for model output -----#
    cdef public float z_delta #height of the energy line - changes for every calculation
    cdef public float flux    #routing flux - might be interpreted as pseudo-mass - between 0 and 1
    cdef public bint isStartCell #boolean variable reflecting if a cell is a release area/start cell
    
    cdef public float min_distance, max_distance #min and max projected path lenght from start cell to current cell
    cdef public float min_gamma, max_gamma #min and max local angle of reach along min_dist path (max_gamma) or max_dist_path (min_gamma)
    cdef public double sl_gamma #straight line travel angle between current cell and start-cell ("geometrical travel angle")
    
    cdef float[:, :] z_delta_neighbour
    cdef float[:, :] z_gamma
    cdef float[:, :] z_alpha
    cdef float[:, :] beta
    cdef float tan_alpha
    
    #-----3x3 arrays characterizing the direct neighbourhood around a cell-----#
    cdef public float[:, :] dem_ng #array/memview with elevation values of the 3x3 neighbourhood
    cdef public float[:, :] tan_beta #local slopes from central cell to neighobur cells in 3x3 neighbourhood
    cdef public float[:, :] dist #horizontal distances from central cell to 3x3 neighbourhood
    cdef public float[:, :] persistence #3x3 persistence weights 
    cdef public float[:, :] distance #3x3 persistence weights
    cdef public short[:, :] no_flow #3x3 neighbourhood of 0|1 masking eligible neighbours -- used to prevent spreading to parent cell
    cdef public float[:, :] r_t #terrain based routing components in 3x3 neighbourhood
    #cdef double[:, :] z_delta_neighbour
    cdef float[:, :] ds
   
    
    cdef public Cell parent #not sure if this is a list or what this should be in flow_class.py --> needs to be checked!!
    cdef public list lOfParents
    cdef public Cell startcell
    cdef public bint is_start
    
    
    def __init__(self,long rowindex, long colindex, float cellsize, 
                 float[:, :] dem_ng,float flux, float z_delta,
                 float alpha, float exp, float flux_threshold, float max_z_delta,
                 bint is_start = False, Cell parent = None, Cell startcell = None):
        '''
        initiation of tan_beta, dist, persistence, no_flow is done on every creation of a new Cell-object in the
        original flow_class code with np.ones_like(3x3) or np.zeros_like(3x3)
        since these 3x3 arrays principally stay the same for every call of __init__ they might be created once
        in flow_core.py and then passed as arguments to this __init__ function - dist can also be calculated once
        outside of flow_class already, since it stays the same for a fixed width grid/raster
        
        e.g.: in flow_core --> 
        
        tan_beta = np.zeros_like(3x3)
        dist = np.array([sqrt(2)*cellsize,  cellsize,   sqrt(2)*cellsize],
                        [cellsize,              0,      cellsize],
                        [sqrt(2)*cellsize,  cellsize,   sqrt(2)*cellsize])
        
        myNewCell = Cell(...,...,dist, tan_beta,...)
        
        NB:
        
        '''
        self.is_start = is_start
        self.rowindex = rowindex
        self.colindex = colindex
        self.cellsize = cellsize
        
        self.dem_ng = dem_ng
        #print(self.dem_ng[0,0])
        self.altitude = self.dem_ng[1, 1]
        
        self.alpha = alpha
        self.exp = exp
        self.flux_threshold = flux_threshold
        self.max_z_delta = max_z_delta
        self.tan_alpha = np.tan(np.deg2rad(self.alpha))
                
        self.tan_beta = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],dtype=np.float32)
        self.dist = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],dtype=np.float32)
        self.persistence = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],dtype=np.float32)
        self.no_flow = np.array([[1., 1., 1.], [1., 1.,1.], [1., 1., 1.]],dtype=np.short)
        self.r_t = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],dtype=np.float32)
        self.ds = np.array([[SQRT2, 1, SQRT2], [1, 1, 1], [SQRT2, 1, SQRT2]],dtype=np.float32)
        self.distance = np.array([[SQRT2*self.cellsize, self.cellsize, SQRT2*self.cellsize], 
                                  [self.cellsize, 0, self.cellsize], 
                                  [SQRT2*self.cellsize, self.cellsize, SQRT2*self.cellsize]],dtype=np.float32)
        
        arrzDN = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],dtype=np.float32) #np.zeros_like(self.tan_beta,dtype=np.float64)
        self.z_delta_neighbour=arrzDN #initializing z_delta_neighbour with zeros
        self.z_gamma=arrzDN
        self.z_alpha=arrzDN
        self.beta=arrzDN
        
        self.flux = flux
        self.z_delta = z_delta
        
        self.FSI = 0.
        self.min_afForest, self.max_afForest = 2.,10.
        self.min_adForest, self.max_adForest = 3e-4, 1e-5
        
        self.maxEffectVForest_F, self.maxEffectVForest_D = 30., 30.
        
        #self.RAD90 = np.deg2rad(90)
        
        '''
        original code allows passing either a variable of type bool or type
        Cell to startcell - this is ambiguous, so in this version is_start(bool) or
        startcell (Cell) are passed as args to the constructor
        '''
        
        if self.is_start == True:
            pass
        else:
            if type(startcell)==Cell:
                self.startcell=startcell
        
        self.lOfParents = []
        if type(parent) == Cell:
            self.lOfParents.append(parent)
            
        
        
        #self.parent=None
        
    def __repr__(self):
        cellStr=f"""
Cell ({self.rowindex}, {self.colindex})
     alpha: {self.alpha}, exp: {self.exp}
     fluxTh:{self.flux_threshold} , maxZDelta: {self.max_z_delta}
     
     flux: {self.flux}
     maxGamma: {self.max_gamma} , minGamma: {self.min_gamma}
     
     
     DEM:
     [[{self.dem_ng[0,0]},{self.dem_ng[0,1]},{self.dem_ng[0,2]}]
     [{self.dem_ng[1,0]},{self.dem_ng[1,1]},{self.dem_ng[1,2]}]
     [{self.dem_ng[2,0]},{self.dem_ng[2,1]},{self.dem_ng[2,2]}]]
     
     Dist:
     [[{self.dist[0,0]},{self.dist[0,1]},{self.dist[0,2]}]
     [{self.dist[1,0]},{self.dist[1,1]},{self.dist[1,2]}]
     [{self.dist[2,0]},{self.dist[2,1]},{self.dist[2,2]}]]
     
     Local Slopes (tan_beta):
     [[{self.tan_beta[0,0]},{self.tan_beta[0,1]},{self.tan_beta[0,2]}]
     [{self.tan_beta[1,0]},{self.tan_beta[1,1]},{self.tan_beta[1,2]}]
     [{self.tan_beta[2,0]},{self.tan_beta[2,1]},{self.tan_beta[2,2]}]]
     
     Persistence:
     [[{self.persistence[0,0]},{self.persistence[0,1]},{self.persistence[0,2]}]
     [{self.persistence[1,0]},{self.persistence[1,1]},{self.persistence[1,2]}]
     [{self.persistence[2,0]},{self.persistence[2,1]},{self.persistence[2,2]}]]
     
     Terrain based routing:
     [[{self.r_t[0,0]},{self.r_t[0,1]},{self.r_t[0,2]}]
     [{self.r_t[1,0]},{self.r_t[1,1]},{self.r_t[1,2]}]
     [{self.r_t[2,0]},{self.r_t[2,1]},{self.r_t[2,2]}]]
     
     No Flow:
     [[{self.no_flow[0,0]},{self.no_flow[0,1]},{self.no_flow[0,2]}]
     [{self.no_flow[1,0]},{self.no_flow[1,1]},{self.no_flow[1,2]}]
     [{self.no_flow[2,0]},{self.no_flow[2,1]},{self.no_flow[2,2]}]]
        """
        return(cellStr)
    
    cpdef add_os(self, double flux):
        '''
        adding flux to a cell
        '''
        self.flux+=flux
        
    cpdef add_parent(self, Cell parent):
        '''
        adding additional parent ot list of parents
        '''
        assert type(parent)==Cell
        self.lOfParents.append(parent)
        
    cpdef calc_fp_travelangle(self):
        '''
        calculates the maximum travel angle to the cell based on the lenght the full path (fp), 
        which is based on minimal self.min_distance of parent cells.
        '''
        cdef list dist_min = []
        cdef double dh, dx, dy
        #dist_min = []
        dh = self.startcell.altitude - self.altitude
        for parent in self.lOfParents:
            dx = cMath.fabs(parent.colindex - self.colindex)
            dy = cMath.fabs(parent.rowindex - self.rowindex)
            dist_min.append(cMath.sqrt(dx ** 2 + dy ** 2) * self.cellsize + parent.min_distance)
        
        if len(dist_min)>0:
            self.min_distance = min(dist_min)
        
        if self.min_distance > 0: #avoiding potential zero division errors
            self.max_gamma = np.rad2deg(np.arctan(dh / self.min_distance))

    cpdef calc_sl_travelangle(self):
        '''
        simple calulation of geometric travel angle from start-cell to current cell
        using shortest distance and altitude diff between cells
        '''
        cdef double dx,dy,dh,ds
        
        dx = cMath.fabs(self.startcell.colindex - self.colindex)
        dy = cMath.fabs(self.startcell.rowindex - self.rowindex)
        dh = self.startcell.altitude - self.altitude

        ds = cMath.sqrt(dx ** 2 + dy ** 2) * self.cellsize
        if ds > 0:
            self.sl_gamma = cMath.atan(dh/ds)*(180./self.PI)
        
        #np.rad2deg(np.arctan(dh / ds))
        
      
    cpdef calc_z_delta(self):
        '''
        funciton calculates the energy line height (z_delta) for all neighbours in the 3x3
        neighbourhood:
        
        self.z_gamma ... ground elevation differences to central cell (values from DEM)
        self.z_alpha ... energy line height differences to central cell (based on alpha and distance to cc) 
        self.z_delta_neighbour ... energy line height for neighbour cell based on the sum of z_delta(current 
                                   elh in cc, the elevation difference (self.z_gamma) and the difference
                                   in energy line height (self.z_alpha)
                                   
        if energy line height is negative (meaning, that the process cannot propagate to neighbour cell based
        on the energy line criterion) the elh in z_delta_neighbour is set to 0, if the z_delta to a neighbour ex
        ceeds the energy line limit (max_z_delta) it is limited to this value
        '''        
        cdef int i,j
        
        for i in range(3):
            for j in range(3):
                #print(self.dem_ng[i,j])
                self.z_gamma[i,j] = self.dem_ng[i,j]-self.altitude
                #print(self.z_gamma[i,j])
                self.z_alpha[i,j] = self.distance[i, j] * self.tan_alpha
                #print(self.z_alpha[i,j])
                #self.ds[i,j] * self.cellsize * self.tan_alpha
                self.z_delta_neighbour[i,j] = self.z_delta + self.z_gamma[i,j] + self.z_alpha[i,j]
                #print('zDNeighbour',i,j,self.z_delta_neighbour[i,j])
                if self.z_delta_neighbour[i,j]<0:
                    self.z_delta_neighbour[i,j] = 0
                if self.z_delta_neighbour[i,j]>self.max_z_delta:
                    self.z_delta_neighbour[i,j]=self.max_z_delta
                #print('zDNeighbour',i,j,self.z_delta_neighbour[i,j])
    
    cpdef calc_tanbeta(self):
        #ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 1, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        #distance = ds * self.cellsize
        cdef int i,j
        cdef double sumTB = 0.
        cdef double sumTB_exp = 0.
                
        for i in range(3):
            for j in range(3):
                
                if (i!= 1) and (j !=1):
                    self.beta[i,j] = cMath.atan((self.altitude-self.dem_ng[i,j])/self.distance[i,j]) + HALFPI
                    self.tan_beta[i,j] = cMath.tan(self.beta[i,j]/2.)
                
                    if self.z_delta_neighbour[i,j] <= 0:
                        self.tan_beta[i,j] = 0
                    if self.persistence[i,j] <= 0:
                        self.tan_beta[i,j] = 0
                else:
                    self.tan_beta[i,j]=0
                
                #print(self.tan_beta[i,j])
                sumTB+=self.tan_beta[i,j]
                sumTB_exp+=self.tan_beta[i,j]**self.exp
                #print(sumTB,sumTB_exp)
        
        #print('#########')
        if cMath.fabs(sumTB) > 0:
            for i in range(3):
                for j in range(3):
                    self.r_t[i,j] = self.tan_beta[i,j]**self.exp/sumTB_exp
                    #print(self.r_t[i,j])
        #print('#########')
            
    cpdef calc_persistence(self):
        cdef int dx,dy
        cdef int idx
        
        #self.persistence = np.zeros_like(self.dem_ng) #no need to define here, can be passed
        if self.is_start==True:
            print('SELF.ISSTART=',self.is_start)
            self.persistence[0,0] += 1
            self.persistence[0,1] += 1
            self.persistence[0,2] += 1
            self.persistence[1,0] += 1
            self.persistence[1,1] += 1
            self.persistence[1,2] += 1
            self.persistence[2,0] += 1
            self.persistence[2,1] += 1
            self.persistence[2,2] += 1
        #elif len(self.lOfParents)>0:
        #    if self.lOfParents[0].is_start==True: #DON'T understand this piece... why only consider parent[0] ??
        #        print('SELF.ISSTART of parents[0]=',self.lOfParents[0].is_start)
        #        self.persistence[0,0] += 1
        #        self.persistence[0,1] += 1
        #        self.persistence[0,2] += 1
        #        self.persistence[1,0] += 1
        #        self.persistence[1,1] += 1
        #        self.persistence[1,2] += 1
        #        self.persistence[2,0] += 1
        #        self.persistence[2,1] += 1
        #        self.persistence[2,2] += 1
        else:
            #print('SELF.ISSTART of parents[0]=',self.lOfParents[0].is_start)
            for idx in range(len(self.lOfParents)):
                #print('parent',idx)
                
            #for parent in self.parent:
                dx = (self.lOfParents[idx].colindex - self.colindex) 
                dy = (self.lOfParents[idx].rowindex - self.rowindex)

                self.no_flow[dy + 1,dx + 1] = 0  # 3x3 Matrix of ones, every parent gets a 0, so no flow to a parent field.
                
                maxweight = self.lOfParents[idx].z_delta
                print('maxweight',maxweight)
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
    
    cpdef calc_distribution(self):
        '''
        calculates distribution of flux to cells of the 3x3 neighbourhood
        '''
        
        cdef int i,j
        cdef double sumRT = 0
        cdef double sumDist = 0
        cdef double sumPxRT = 0
        cdef int count = 0
        cdef int countEligible = 0
        cdef double massToDistribute = 0
        
        self.calc_z_delta()                                               #note: seems to work
        self.calc_persistence()
        for i in range(3):
            for j in range(3):
                self.persistence[i,j]*=self.no_flow[i,j]
                #print(self.persistence[i,j])
        self.calc_tanbeta()
        
        if self.is_start==False:
            self.calc_fp_travelangle()
            self.calc_sl_travelangle()
        
        threshold = self.flux_threshold
            
        for i in range(3):
            for j in range(3):
                sumRT+=self.r_t[i,j]
                sumPxRT+=self.r_t[i,j]*self.persistence[i,j]
        
        #print('----')
        #print(sumRT)
        #print(sumPxRT)
        #print('----')
         
        if sumRT>0:
            for i in range(3):
                for j in range(3):
                    self.dist[i,j]=(self.persistence[i,j] * self.r_t[i,j])/sumPxRT * self.flux
                    #print(self.dist[i,j])
            
        #distribution of flux from cells with flux > 0, but flux < threshold value Rstop
        #count #number of cells which fall below the threshold but are larger than 0
        
        for i in range(3):
            for j in range(3):
                if (self.r_t[i,j]>=threshold):
                    countEligible+=1
                if (self.r_t[i,j]>0) and (self.r_t[i,j]<threshold):
                    count+=1
                    massToDistribute+=self.r_t[i,j]
        
        print(countEligible)
        
        if (massToDistribute>0) and (count>0):
            '''A) wenn noch mass to distribute, dann verteile sie gleichmäßig auf
               andere Nachbarzellen, die >= threshold sind'''
            for i in range(3):
                for j in range(3):
                    if self.dist[i,j]>=threshold:
                        self.dist[i,j]+=massToDistribute/count
                        sumDist+=self.dist[i,j]
                    else:
                        self.dist[i,j]=0
        
        if (sumDist < self.flux) and (count >0):
            '''B) verstehe ich noch nicht ganz,
               was diese funktion macht ...
               wohl sowas wie:
                   -wenn trotzdem noch flux auf der Zelle ist, dann teil die
                    Differenz nochmal gleich auf wie unter A) oben'''
            for i in range(3):
                for j in range(3):
                    if self.dist[i,j]>=threshold:
                        self.dist[i,j]+=(self.flux-sumDist)/count
        
        cdef np.ndarray[long, ndim=1] row_local = np.zeros((countEligible,),dtype=np.int)
        cdef np.ndarray[long, ndim=1] col_local = np.zeros((countEligible,),dtype=np.int)
        
        cdef np.ndarray[long, ndim=1] row_return = np.zeros((countEligible,),dtype=np.int)
        cdef np.ndarray[long, ndim=1] col_return = np.zeros((countEligible,),dtype=np.int)
        
        cdef np.ndarray[float, ndim=1] dist_return = np.zeros((countEligible,),dtype=np.float32)
        cdef np.ndarray[float, ndim=1] zDn_return = np.zeros((countEligible,),dtype=np.float32)
        
        cdef int idx1d = 0
        for i in range(3):
            for j in range(3):
                #print(self.dist[i,j])
                #print(self.dist[i,j],threshold)
                if self.dist[i,j]>=threshold:
                    row_local[idx1d]=i
                    row_local[idx1d]=j
                    
                    row_return[idx1d]=self.rowindex-1+i
                    col_return[idx1d]=self.colindex-1+j
                    
                    dist_return[idx1d]=self.dist[i,j]
                    zDn_return[idx1d] =self.z_delta_neighbour[i,j]
                    idx1d+=1
                    
                
        #row_local, col_local = np.where(self.dist > threshold)
        return row_return,col_return,dist_return,zDn_return
