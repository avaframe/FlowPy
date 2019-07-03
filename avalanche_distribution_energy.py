import numpy as np
import sys
import time

sys.path.append('/home/W/Neuhauser/Frei/python_libs/')
import raster_io as io

class Cell():
    
    def __init__(self, rowindex, colindex, altitude, dem_ng, cellsize, startcell):
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = altitude
        self.dem_ng = dem_ng
        self.cellsize = cellsize
        self.calc_tanbeta()
        self.alpha = 25
        if startcell == True: #check, if start cell exist (start cell is release point)
            self.is_start = True # set is_satrt to True
        else:            
            self.startcell = startcell # set is_satrt to True
            self.is_start = False # set is_satrt to False

        self.calc_energy()
        
    def calc_energy(self):
        if self.is_start:
            self.elh = 0
        else:
            dx = np.abs(self.startcell.colindex - self.colindex) * self.cellsize
            dy = np.abs(self.startcell.rowindex - self.rowindex) * self.cellsize
            ds = np.sqrt(dx**2 + dy**2)
            dz = self.startcell.altitude - self.altitude
            
            self.elh = dz - ds * np.tan(np.deg2rad(self.alpha))
        
    def steepest_descend(self):
        min_alt = np.amin(self.dem_ng)
        row_idx, col_idx = np.where(self.dem_ng == min_alt)
        row_idx = row_idx[0]
        col_idx = col_idx[0]
        out_row = self.rowindex - 1 + row_idx
        out_col = self.colindex - 1 + col_idx
        return out_row, out_col
    
    def calc_tanbeta(self):
        self.tan_beta = np.zeros_like(self.dem_ng)
        for i in range(3):
            for j in range(3):
                if i == 1 or j == 1:
                    distance = self.cellsize
                else:
                    distance = np.sqrt(2) * self.cellsize
                self.tan_beta[i, j] = (self.dem_ng[i, j] - self.altitude)/distance
            self.tan_beta[self.tan_beta > 0] = 0
            
    def calc_distribution(self):
        self.dist = np.zeros_like(self.dem_ng)
        exp = 8
        threshold = 0.01
        for i in range(3):
            for j in range(3):
                self.dist[i, j] = self.tan_beta[i, j] ** exp / np.sum(self.tan_beta ** exp)
        row_local, col_local = np.where(self.dist > threshold)
        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local
        
#Reading in the arrays
path = '/home/P/Projekte/18130-GreenRisk4Alps/Simulation/PAR6_ValsGries_AUT/'
file = path + 'dhm_cropped.asc'
release_file = path + 'init/init_ava_c1.asc'
ava_out = path + 'distribution_gires.asc'
elh_out = path + 'energy_gries.asc'

dem, header = io.read_raster(file)
cellsize = header["cellsize"]   
release, header = io.read_raster(release_file) 
ava_steep = np.zeros_like(dem)
elh = np.zeros_like(dem)

start = time.time()


cell_list = []  
row_list, col_list = np.where(release == 1)

a = list(range(0, len(row_list)))

for i in a:
    print("Calculating Startcell: " + str(i) + " from " + str(a[-1]))
    cell_list = []
    row_idx = row_list[i]
    col_idx = col_list[i]
    
    ava_steep[row_idx, col_idx] = 1
    
    dem_ng = dem[row_idx - 1 : row_idx + 2, col_idx - 1 : col_idx + 2]
    if np.size(dem_ng) < 9:
        continue
    
    startcell = Cell(row_idx, col_idx, dem[row_idx, col_idx], dem_ng, cellsize, startcell=True)
    
    cell_list.append(startcell)
    checked = 0
    for cells in cell_list:
        if cells.elh < 0: #Stopping Condition
            continue
        
        row, col = cells.calc_distribution()
        for i in range(checked, len(cell_list)):
            k = 0
            while k < len(row):
                if row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex:
                    row = np.delete(row, k)
                    col = np.delete(col, k)
                else:                                 
                    k += 1
                    
        for k in range(len(row)):
            dem_ng = dem[row[k]-1:row[k]+2,col[k]-1:col[k]+2]
            if np.size(dem_ng) < 9:
                continue
            cell_list.append(Cell(row[k], col[k], dem[row[k], col[k]], dem_ng, cellsize, cell_list[0]))
            ava_steep[row[k], col[k]] += 1
        checked += 1 
        
        elh[cells.rowindex, cells.colindex] = max(elh[cells.rowindex, cells.colindex], cells.elh)
        
                
end = time.time()            
print('Time needed: ' + str(end - start) + ' seconds')

io.output_raster(file, ava_out, ava_steep)
io.output_raster(file, elh_out, elh)
