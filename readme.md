# Readme

## Input Files

All input files need to be .asc or .tif files. 
All files need the same raster resolution, normal sizes are 5x5 or 10x10 meters. 
All Layers need the exact same extend. If not the Code will give you feedback which layer is not acurate.

### Input Files:

- DEM:
	- The DEM Raster file needs a no data value lower then zero. standard = -9999
	
- Release Zones:
	- The release layer needs values higher then zero for the release pixels. 
	
### Optional Input Files:

- Infrastructure:
	- The inrastructure Layer needs values higher then zero. Different values can be used for different infrastructures.
	- The backcalculation layer has the information of the infrastructure that was hit. Higher values win over lower ones.
	
## Output

- Energy:

- Cell Counts:

## Code

./main.py 25 8 Avalanche ./examples/helix/ ./examples/helix/dhm.asc ./examples/helix/class_1.asc
