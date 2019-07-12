import numpy as np
import sys
import time

sys.path.append('/home/W/Neuhauser/Frei/python_libs/')
import raster_io as io

class Cell():
    
    def __init__(self, rowindex, colindex, dem_ng, cellsize, mass, startcell):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = dem_ng[1,1]
        self.dem_ng = dem_ng
        self.cellsize = cellsize
        #self.tan_beta = np.zeros_like(self.dem_ng)
        #self.dist = np.zeros_like(self.dem_ng)
        self.alpha = 36
        #self.velocity_sqr = 0
        self.mass = mass
        self.elh_ng = np.ones_like(dem_ng)
        if startcell == True: #check, if start cell exist (start cell is release point)
            self.is_start = True # set is_satrt to True
        else:            
            self.startcell = startcell # set is_satrt to True
            self.is_start = False # set is_satrt to False        
        #self.calc_tanbeta()
        self.calc_energy()
        self.calc_energy_neighbour()
        self.calc_distr_projection()
        #self.calc_mass()
        
    def calc_energy(self):
        if self.is_start:
            self.elh = 0
        else:
            dx = np.abs(self.startcell.colindex - self.colindex) * self.cellsize #g
            dy = np.abs(self.startcell.rowindex - self.rowindex) * self.cellsize
            ds = np.sqrt(dx**2 + dy**2)
            dz = self.startcell.altitude - self.altitude            
            self.elh = dz - ds * np.tan(np.deg2rad(self.alpha))
            
    def calc_energy_neighbour(self):
        #ToDo: check if elh of cell is near 0, otherwise dont do this crap
        if self.is_start:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    dz = self.altitude - dem_ng[i+1, j+1]
                    dx = np.abs(j) * self.cellsize
                    dy = np.abs(i) * self.cellsize
                    ds = np.sqrt(dx**2 + dy**2)
                    self.elh_ng[i+1, j+1] = dz - ds * np.tan(np.deg2rad(self.alpha))
            self.elh_ng[self.elh_ng < 0 ] = 0
            self.elh_ng[self.elh_ng > 0 ] = 1
        else:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    dz = self.startcell.altitude - dem_ng[i+1, j+1]
                    dx = np.abs(self.startcell.colindex - self.colindex + j) * self.cellsize
                    dy = np.abs(self.startcell.rowindex - self.rowindex + i) * self.cellsize
                    ds = np.sqrt(dx**2 + dy**2)
                    self.elh_ng[i+1, j+1] = dz - ds * np.tan(np.deg2rad(self.alpha))
            self.elh_ng[self.elh_ng < 0 ] = 0
            self.elh_ng[self.elh_ng > 0 ] = 1
            
    def calc_velocity(self):# g
        if self.is_start:
            self.elh = 0
        else:
            dx = (self.startcell.colindex - self.colindex) * self.cellsize #g
            dy = (self.startcell.rowindex - self.rowindex) * self.cellsize
            ds = np.sqrt(dx**2 + dy**2)
            dz = self.startcell.altitude - self.altitude            
            self.velocity_sqr = 2 * 9.81 * (dz - ds * np.tan(np.deg2rad(self.alpha)))
            
    def add_mass(self, mass):
        #print(self.mass, mass)
        self.mass = self.mass + mass
        #print('Added mass: ',self.mass)
        
    def steepest_descend(self):
        min_alt = np.amin(self.dem_ng)
        row_idx, col_idx = np.where(self.dem_ng == min_alt)
        row_idx = row_idx[0]
        col_idx = col_idx[0]
        out_row = self.rowindex - 1 + row_idx
        out_col = self.colindex - 1 + col_idx
        return out_row, out_col
    
    def calc_tanbeta(self):        
        for i in range(3):
            for j in range(3):
                if i == 1 or j == 1:
                    distance = self.cellsize
                else:
                    distance = np.sqrt(2) * self.cellsize
                self.tan_beta[i, j] = (self.dem_ng[i, j] - self.altitude)/distance
            self.tan_beta[self.tan_beta > 0] = 0
            
    def calc_distribution(self):
        exp = 8
        threshold = 0.05
        for i in range(3):
            for j in range(3):
                self.dist[i, j] = self.tan_beta[i, j] ** exp / np.sum(self.tan_beta ** exp)
        row_local, col_local = np.where(self.dist > threshold)
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local
    
    def calc_distr_projection(self):  # global
        
        if self.is_start:
            #slopes = self.dem_ng - self.altitude
            # mass stuff
            self.neighbour_cells = self.dem_ng - self.altitude  # because elh=0 neighbour_cells = slopes
            #normalization = np.sum(self.neighbour_cells[self.neighbour_cells < 0])
            #self.neighbour_cells = self.neighbour_cells/normalization
            #row_local, col_local = np.where(self.neighbour_cells > 0)
            #return self.rowindex - 1 + row_local, self.colindex - 1 + col_local
            
        else: 
            #global stuff
            dx = (self.colindex - self.startcell.colindex) * self.cellsize # x component of avalanche flow direction, global
            dy = (self.rowindex - self.startcell.rowindex) * self.cellsize # y component of avalanche flow direction, global
            avi_direction = np.rad2deg(np.arctan2(dy, dx))  # avalanche direction in degrees, global
            #neighbour_direction = np.linspace(0.0, 315, 8)  # direction in deg to all cells
            
            # Setting the direction for the Center Cell in the opositre direction for the avalanche, so it´s not calculated
            neighbour_direction = np.array([[135, 90, 45], [180, np.nan, 0], [225, 270 , 315]])  # direction in deg to all cells, local
            delta_deg = neighbour_direction - avi_direction  # difference between ava direction and the neighbor cells
            direction_projection = np.cos(np.deg2rad(delta_deg)) # direction_projection[1] = NE, direction_projection[2] = N ..., global influence
            
            #global and local
            #neighbour_cells = -self.elh * direction_projection + terrain_slopes[:, cC[0], cC[1]]
            slopes = self.dem_ng - self.altitude
            self.neighbour_cells = self.elh * direction_projection + slopes  # slopes is local influence, add negative scaled projection
            
            #row_local, col_local = np.where(self.neighbour_cells > 0) 
            #return self.rowindex - 1 + row_local, self.colindex - 1 + col_local
        
    def calc_mass(self):
        # mass stuff
        dummy = self.neighbour_cells
        self.neighbour_cells *= self.elh_ng
        normalization = np.sum(self.neighbour_cells[self.neighbour_cells < 0])
        self.neighbour_cells = self.neighbour_cells/normalization# Normilzing the negative values and mutlipling it with the mass of the cell, so they are positive
        dummy2 = self.neighbour_cells
        # print(three_flow)
        #Fdir = mapping(three_flow)
        row_local, col_local = np.where(self.neighbour_cells > 0) 
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.neighbour_cells[row_local, col_local] * self.mass
    
def get_start_idx(release):
    row_list, col_list = np.where(release > 0)  # Gives back the indices of the release areas
    altitude_list = []
    for i in range(len(row_list)):
        altitude_list.append(dem[row_list[i], col_list[i]])    
    altitude_list, row_list, col_list = list(zip(*sorted(zip(altitude_list, row_list, col_list), reverse=True)))  #Sort this lists by altitude
    return row_list, col_list    
        
#Reading in the arrays
path = '/home/neuhauser/git_rep/graviclass/'
file = path + 'dhm.asc'
release_file = path + 'class_1.asc'
elh_out = path + 'energy_newtry.asc'
mass_out = path + 'mass_newtry.asc'

dem, header = io.read_raster(file)
cellsize = header["cellsize"]   
release, header = io.read_raster(release_file) 
elh = np.zeros_like(dem)
mass_array = np.zeros_like(dem)

start = time.time()

cell_list = []  
row_list, col_list = get_start_idx(release)

startcell_idx = 0
while startcell_idx < len(row_list):
    print("Calculating Startcell: " + str(startcell_idx) + " of " + str(len(row_list)) + " = " + str(round(startcell_idx/len(row_list)*100,2)) + "%")
    cell_list = []
    row_idx = row_list[startcell_idx]
    col_idx = col_list[startcell_idx]    
    dem_ng = dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2] # neighbourhood DEM
    if np.size(dem_ng) < 9:
        startcell_idx += 1
        continue
    
    startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1, startcell=True)
    # If this is a startcell just give a Bool to startcell otherwise the object startcell
    
    cell_list.append(startcell)
    checked = 0
    for cells in cell_list:
        if cells.elh < 0: # Stopping Criteria, or cells.mass < 1e-5
            #checked += 1
            continue
        #start_dist =  time.time()
        #row, col = cells.calc_distribution()
        #row, col = cells.calc_distr_projection()
        row, col, mass = cells.calc_mass()
        #end_dist = time.time()
        #print("Distribution took: " + str(end_dist - start_dist) + "seconds" )
        #start_check = time.time()
        for i in range(len(cell_list)):  # Check if Cell already exists
            k = 0
            while k < len(row):
                if (row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex):
                    cell_list[i].add_mass(mass[k])
                    row = np.delete(row, k)
                    col = np.delete(col, k)
                    #print('from cell',cells.rowindex, cells.colindex, 'to cell' , str(cell_list[i].rowindex) , str(cell_list[i].colindex))                    
                    mass = np.delete(mass, k)
                else:                                 
                    k += 1
        #end_check = time.time()
        #print("Check took: " + str(end_check - start_check) + "seconds" )           
        for k in range(len(row)):
            dem_ng = dem[row[k]-1:row[k]+2, col[k]-1:col[k]+2]  # neighbourhood DEM
            if np.size(dem_ng) < 9:
                #checked += 1# Dirty way to don´t care about the edge of the DEM
                continue
            cell_list.append(Cell(row[k], col[k], dem_ng, cellsize, mass[k], startcell))            
        #checked += 1         
        elh[cells.rowindex, cells.colindex] = max(elh[cells.rowindex, cells.colindex], cells.elh)
        mass_array[cells.rowindex, cells.colindex] = cells.mass
        #cell_list.pop(0)
    release[elh > 0] = 0  # Check if i hited a release Cell, if so set it to zero and get again the indexes of release cells
    #for cell in cell_list:
        #mass_array[cell.rowindex, cell.colindex] = cell.mass
    row_list, col_list = get_start_idx(release)
    startcell_idx += 1
        
end = time.time()            
print('Time needed: ' + str(end - start) + ' seconds')

io.output_raster(file, elh_out, elh, 4326)
io.output_raster(file, mass_out, mass_array, 4326)
