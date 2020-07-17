# Readme

Right now the Code only works on Linux Machines, due to the multiprocessing.\
Terminal and Windows Version are in the making...
Run the Code via the main.py script. \
Some PyQt libraries are needed and rasterio.\
There is a .yml file in the repo that includes all needed libraries, you can import this file in anaconda.\
Right now only GUI version is working, working on terminal version.

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