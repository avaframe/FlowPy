import numpy as np
import sys
import time

sys.path.append('/home/W/Neuhauser/Frei/python_libs/')
import raster_io as io

class Cell():
    
    def __init__(self, rowindex, colindex, dem_ng, cellsize, mass, kin_e, parent, startcell):
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
        self.exp = 4 # Exp should be a odd number, otherwise in tan_beta negative values will get positive
        self.mass_threshold = 3*10**-4
        #self.velocity_sqr = 0
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
        self.calc_tanbeta()

        
    def calc_kinetic_energy(self):
        delta_e_kin_pot = (self.dem_ng - dem_ng[1,1]) * (-1)
        ds = np.array([[np.sqrt(2),1,np.sqrt(2)],[1,0,1],[np.sqrt(2),1,np.sqrt(2)]])
        e_friction = ds * self.cellsize * np.tan(np.deg2rad(self.alpha))
        self.kin_energy_neighbour = self.kin_e + delta_e_kin_pot - e_friction
        self.kin_energy_neighbour[self.kin_energy_neighbour < 0] = 0
                    
    def add_mass(self, mass):
        self.mass += mass
        
    def add_parent(self, parent):
        self.parent.append(parent)
           
    def calc_tanbeta(self):  
        exp = self.exp
        dh = self.kin_e/(9.81*self.mass*self.cellsize**2 * 2 * 100)
        print(dh)
        ds = np.array([[np.sqrt(2),1,np.sqrt(2)],[1,0,1],[np.sqrt(2),1,np.sqrt(2)]])
        distance = ds * self.cellsize
        #self.tan_beta = ((self.dem_ng - self.altitude) * (-1)) / distance
        for i in range(3):
            for j in range(3):
                #if ((self.dem_ng[i, j] - (self.altitude + dh)) * (-1)) > 0:
                self.tan_beta[i, j] = (self.dem_ng[i, j] - (self.altitude + dh)) * (-1)/distance[i, j]  # Downslope = postive value
                #else:
                #    self.tan_beta[i, j] = 0.0637 * np.arctan(distance/((self.dem_ng[i, j] - (self.altitude + dh)) * (-1))) + 0.1
            self.tan_beta[self.tan_beta < 0] = 0
        self.tan_beta[self.kin_energy_neighbour <= 0] = 0
        self.tan_beta[self.direction <= 0] = 0
        self.tan_beta[1,1] = 0
        #self.tan_beta = np.abs(self.tan_beta)
        if np.sum(self.tan_beta > 0):
            self.p_fd = self.tan_beta ** exp/ np.sum(self.tan_beta ** exp)
        
    def calc_direction(self):
        if self.is_start:
            self.direction += 1
        elif self.parent[0].is_start:
            self.direction += 1
        else:
            #kin_e_sum = 0
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

                #kin_e_sum += parent.kin_e
            np.rot90(self.direction,2)
            #self.direction /= kin_e_sum
            
                    
    def calc_distribution(self):
        #exp = self.exp
        threshold = self.mass_threshold
        if np.sum(self.p_fd > 0):
            for i in range(3):
                for j in range(3):
                    self.dist[i, j] = self.direction[i, j] * self.p_fd[i, j] / np.sum(self.direction * self.p_fd) * self.mass
                #self.dist = (self.direction * self.kin_energy_neighbour)/np.sum(self.direction * self.kin_energy_neighbour) * self.mass
                #print(self.direction[i, j] * self.tan_beta[i, j] ** exp)
                #print(np.sum(self.direction * self.tan_beta ** exp) * self.mass)
                #print(self.direction[i, j] * self.tan_beta[i, j] ** exp/np.sum(self.direction * self.tan_beta ** exp) * self.mass)
                
        #count = len(self.dist[0 < self.dist < threshold])
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
    
## Programm:        
def get_start_idx(release):
    row_list, col_list = np.where(release > 0)  # Gives back the indices of the release areas
    if len(row_list) > 0:
        altitude_list = []
        for i in range(len(row_list)):
            altitude_list.append(dem[row_list[i], col_list[i]])    
        altitude_list, row_list, col_list = list(zip(*sorted(zip(altitude_list, row_list, col_list), reverse=True)))  #Sort this lists by altitude
    return row_list, col_list   

# =============================================================================
# def back_calculation(cell):
#     back_list = []
#     parent = cell.parent[0]
#     while parent:
#         for parent in cell.parent:
#             back_list.append(parent)
# =============================================================================
        
#Reading in the arrays
path = '/home/neuhauser/git_rep/graviclass/'
file = path + 'dhm.asc'
release_file = path + 'class_1.asc'
# =============================================================================
# path = '/home/P/Projekte/18130-GreenRisk4Alps/Simulation/PAR3_Oberammergau/'
# file = path + 'DEM_10_3.tif'
# release_file = path + 'init/release_class_1.asc'
# =============================================================================
elh_out = path + 'energy_flowr.asc'
mass_out = path + 'mass_flowr.asc'
index_out = path + 'index_flowr.asc'

#dem = np.array([[920,915,920], [910,900,910], [880, 890, 885]])
#cellsize = 10
#release = np.array([[0,0,0], [0,1,0], [0,0,0]])

dem, header = io.read_raster(file)
cellsize = header["cellsize"]  
release, header = io.read_raster(release_file) 

elh = np.zeros_like(dem)
mass_array = np.zeros_like(dem)
index_array = np.zeros_like(dem)


start = time.time()
#Core
cell_list = []  
row_list, col_list = get_start_idx(release)

startcell_idx = 0
while startcell_idx < len(row_list):
    sys.stdout.write('\r' "Calculating Startcell: " + str(startcell_idx+1) + " of " + str(len(row_list)) + " = " + str(round(startcell_idx+1/len(row_list)*100,2)) + "%" '\r')
    sys.stdout.flush()
    #print("Calculating Startcell: " + str(startcell_idx+1) + " of " + str(len(row_list)) + " = " + str(round(startcell_idx/len(row_list)*100,2)) + "%")
    cell_list = []
    row_idx = row_list[startcell_idx]
    col_idx = col_list[startcell_idx]    
    dem_ng = dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2] # neighbourhood DEM
    if np.size(dem_ng) < 9:
        startcell_idx += 1
        continue
    
    startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1, 15, None, startcell=True)
    # If this is a startcell just give a Bool to startcell otherwise the object startcell
    
    cell_list.append(startcell)
    checked = 0
    index = 0
    for cells in cell_list:
        row, col, mass, kin_e = cells.calc_distribution()

        for i in range(int(checked), len(cell_list)):  # Check if Cell already exists
            k = 0
            while k < len(row):
                if (row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex):
                    cell_list[i].add_mass(mass[k])
                    cell_list[i].add_parent(cells)
                    cell_list[i].calc_direction()
                    row = np.delete(row, k)
                    col = np.delete(col, k)
                    mass = np.delete(mass, k)
                    kin_e = np.delete(kin_e, k)
                else:                                 
                    k += 1
         
        for k in range(len(row)):
            dem_ng = dem[row[k]-1:row[k]+2, col[k]-1:col[k]+2]  # neighbourhood DEM
            if np.size(dem_ng) < 9:
                #checked += 1# Dirty way to don´t care about the edge of the DEM
                continue
            cell_list.append(Cell(row[k], col[k], dem_ng, cellsize, mass[k], kin_e[k], cells, startcell))            
        #checked += 1         
        elh[cells.rowindex, cells.colindex] = max(elh[cells.rowindex, cells.colindex], cells.kin_e)
        mass_array[cells.rowindex, cells.colindex] = cells.mass
        index_array[cells.rowindex, cells.colindex] = index
        index += 1
    release[elh > 0] = 0  # Check if i hited a release Cell, if so set it to zero and get again the indexes of release cells

    row_list, col_list = get_start_idx(release)
    startcell_idx += 1
        
end = time.time()            
print('Time needed: ' + str(end - start) + ' seconds')

# Output
io.output_raster(file, elh_out, elh, 4326)
io.output_raster(file, mass_out, mass_array, 4326)
io.output_raster(file, index_out, index_array, 4326)
