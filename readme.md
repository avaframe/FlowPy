#Flow Py

To define protective functions and to quantify the protective effects of forests, we created the Flow-Py model that 
identifies process areas of gravitational hazards, including avalanches, rockfall and debris slides. The model is 
written in Python to keep it easy adjustable. The run out routine of Flow-Py is based on the principles of energy 
conservation including frictional dissipation assuming simple coulomb friction, leading to constant travel-angle. 
Potential release areas and the corresponding travel angle have to be adapted for each type of gravitational mass movements. 
A important improvement, compared to similar models, is that it can handle mass movement in flat and uphill terrain. 
One major advantage of this model is its simplicity, resulting in a computationally inexpensive implementation, which 
allows for an application on a regional scale, covering large simulation areas. The adaptivity of the model further 
allows to consider existing infrastructure and to detect starting zones endangering the corresponding areas in a back-calculation step. 
Additionally, by adding forest cover to the simulations we can identify which forest area has a protective function and, 
based on information about forest structure, calculate the protective effect this forest provides to down slope infrastructure.

## Running the Code
./main.py --gui -> Gui Version\
./main.py alpha_angle exponent process working_directory path_to_dem path_to_release path_to_infrastructure \
- alpha_angle: max. runout angle for the process: Austria -> Avalanche 25, Rockfall 35, Debris Slides 22 \
- exponent: controls the lateral spreading, avalanches 8, rockfall and debris slides 75 (= single flow, except in flat terrain)\
- process: can be -> Avalanche Rockfall or Soil Slides, controls max. energy line height\
- working_directory: where to create and save result folder\
- path_to_dem: well it's the path to the DEM (.asc or .tif)
- path_to_release: well it's the path to the release layer\
- path_to_infrastructure: well it's the path to the Infra layer, OPTIONAL!!!\

Right now the Code only works on Linux Machines, due to the multiprocessing.\
Run the Code via the main.py script.\
Some PyQt libraries are needed and rasterio.\
There is a .yml file in the repo that includes all needed libraries, you can import the environment in anaconda.\

## Input Files

All input files need to be .asc or .tif files. \
All files need the same raster resolution, normal sizes are 5x5 or 10x10 meters. \
All Layers need the exact same extend. If not, the Code will give you feedback which layer is not accurate.

### Input Files:

- DEM:
	- The DEM Raster file needs a no data value lower then zero. standard = -9999
	
- Release Zones:
	- The release layer needs values higher then zero for the release pixels. NoData < 0, or -9999
	
### Optional Input Files:

- Infrastructure:
	- The infrastructure layer needs values higher then zero for infrastructure. Different values can be used for 
	different infrastructures classes and will be saved in the backcalculation.
	- The backcalculation layer has the information of the infrastructure that was hit. Higher values win over lower ones.
	
## Output

- Energy line height:
    Includes the highest energy line height for every pixel.
- Sum of energy line height:
    The energy line height summed up on every pixel.
- Susceptibility:
    The result of the susceptibility calculation for every pixel.
- Cell Counts:
    Saves how often one pixel gets a hit from different release points.
- Backcalculation:
    Saves the backcalculation from Infrastructure to release point

## Contact

For Questions contact: \
Michael Neuhauser, Austrian Research Centre for Forest: Michael.Neuhauser@bfw.gv.at \
Christopher D'Amboise, Austrian Research Centre for Forest: Christopher.DAmboise@bfw.gv.at