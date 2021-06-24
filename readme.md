# Flow-Py

Flow-Py is an adaptable open source implementation of existing approaches to GMM (gravitational mass movements) run out models. The main objective of this tool is to compute the spatial run out distribution (routing and stopping) for GMMs in three dimensional terrain. This model has been designed to be computationally light, allowing it application on a regional scale. 

The data based (empirical) modeling approach of Flow-Py, combining routing and stopping run out in three dimensional terrain, allows to identify process areas and corresponding magnitudes. The model is motivated by a combination of existing process and data based approaches, requiring a minimum of input data, providing valuable output. By considering the spatial evolution (no temporal evolution equations are solved in the model) the resulting run out is mainly dependent by the terrain and the location of the starting point

The model is written in Python to keep it easy adjustable. The run out routine of Flow-Py is based on the principles of energy conservation including frictional dissipation assuming simple coulomb friction, leading to constant travel-angle. Potential release areas and the corresponding travel angle have to be adapted for each type of gravitational mass movements. 
A important improvement, compared to similar models, is that it can handle mass movement in flat and uphill terrain. One major advantage of this model is its simplicity, resulting in a computationally inexpensive implementation, which allows for an application on a regional scale, covering large simulation areas. The adaptivity of the model further allows to consider existing infrastructure and to detect starting zones endangering the corresponding areas in a back-calculation step. 

## Running the Code

python3 main.py --gui -> Gui Version  
python3 main.py alpha_angle exponent working_directory path_to_dem path_to_release infra=path_to_infrastructure(Optional) flux_threshold=positiv_number(Optional) max_z_delta=positiv_number(Optional)

- alpha_angle: max. runout angle for the process: Austria -> Avalanche 25, Rockfall 35, Debris Slides 22 
- exponent: controls the lateral spreading, avalanches 8, rockfall and debris slides 75 (= single flow, except in flat terrain) , an exponent of 1 delivers wide spreading while when we increase it towards infinity we get a single flow
- working_directory: where to create and save result folder 
- path_to_dem: well it's the path to the DEM (.asc or .tif)
- path_to_release: well it's the path to the release layer 
- path_to_infrastructure: well it's the path to the Infra layer, Optional!
- flux_threshold: when Flow-Py stops the spreading, Standard = 0.0003, Optional!
- max_z_delta: The max. z_delta your process can reach. Some hints: Avalanche = 270 /  Rockfall = 50 / Soil Slides = 12 / Standard = 8848 (no limitation), Optional!

Right now the code works on Linux and Windows machines. Haven't tested it on MacOS, if you are able to run it there, please give us feedback.
Run the Code via the main.py script: python3 main.py ...  
Some PyQt libraries are needed and rasterio.  
There is the requirements.txt file in the repo that includes all needed libraries. (Work in progress...)

If you have trouble with GDAL or rasterio on Windows use this links to get the needed version directly from their website, first install GDAL and then rasterio.

GDAL: https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal

rasterio: https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio

## Input Files

All input raster files are required to be .asc or .tif files.  
All files need the same raster resolution, normal sizes are 5x5 or 10x10 meters.  
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

- z_delta:
    Includes the highest z_delta for every cell.
- Sum of z_delta:
    z_delta summed up on every cell.
- Flux:
    The result of the flux calculation for every cell.
- Cell Counts:
    Saves how often one pixel gets a hit from different release points.
- Backtracking:
    Saves the backtracking from Infrastructure to release point
- Flow Path Travel Angle, FP_TA:
    Saves the &gamma; angle along the flow path
- Straight Line Travel Angle, SL_TA:
    Saves the &gamma; angle, while the distances are calculated via a straight line from the release cell to the current cell

## Motivation 2D



![Image](img/Motivation_2d.png)
*Fig. 1: Definition of angles and distances for the calculation of z_delta, where s is the distance along the path and z(s) the corresponding altitude.*

The model equations that control the run out in three dimensional terrain are mostly motivated with respect to geometric two dimensional concepts, that control the main routing and final stopping of the flow.

Figure 1 summarizes the basic concept of a constant run out angle (alpha) with the corresponding geometric relations in two dimensions along a possible process path.

![tan_alpha](img/tan_alpha.png)

The local travel angle gamma is defined by the height and distance difference from the release point to the current location.

![tan_gamma](img/tan_gamma.png)

The angle delta is the difference between gamma and alpha and is associated to the energy left in the process, so when delta equals zero or gamma equals alpha, the maximum runout distance is reached.

![z_alpha](img/z_alpha.png)

Z_alpha can be interpreted as dissipation energy.

![z_gamma](img/z_gamma.png)

Z_gamma is the height difference between the starting point and the current calculation step at s. 

Z_delta is the difference between Z_gamma and Z_alpha, so when Z_delta is lower or equal zero the calculation stops. Z_delta can also be interpreted as the energy left in the process.

![z_delta](img/z_delta.png)

## Handling the Spatial Input

To run the model at least 2 main raster inputs are needed. First the digital elevation model on which we solve the equations, and second the release raster, which defines were the starting points or release cells are on the raster. 

For every release cell we start the calculation. The goal is to determine potential children via 2 stopping criteria. 

- Z_delta has to be greater then zero: Z_delta > 0
- Routing Flux has to be greater then the Routing cut off: R_i > R_Stop

If the criteria are fulfilled the child is saved in the path and for every cell/child in the path we search for new potential children.

If we reach the end of the path and no new children fulfill the criteria we close the calculation and save the needed information from the path to our output raster files. Then the calculation starts again for the next release cell. The spatial extend and magnitude for all release cells are summarized in the output raster files, which represent the overlay of all paths.

Every path is independent from the other, but depending on the information we want to extract, we save the highest values (Z_delta) or sums (Cell Counts) of different paths to the output raster file.

## Calculation Steps on the Path

Here we will go through all calculation steps as they are in the code under: flow_class.calc_distribution()

To bring the thoughts from the motivation from 2D model to a 3D grid we must implement a few new definitions.

​	![grid_overview](img/Neighbours.png)

First we need to bring in the definition for base. This is the current raster cell we are looking at, and from which we do our calculations. This would be at distance s along the path in Fig. 1. 

Every base can have one or more parents, except the starting cell, where we start our calculation, this would be at s = s_0 in Fig. 1.

### z_delta

From the base we solve now the equations (6,7 and 8) for every neighbor n, if Z_bn^{delta} is higher then zero, this neighbor is defined as a potential child of this base, and spreading is allowed in this direction.

![z_delta_i](img/z_delta_array.png)

Here S_bn is the projected distance between the base and the neighbor.

As Z_bn^{delta} can be interpreted as  kinetic energy it is possible to limit this value to a maximum. Regarding physical models this would correspond to a turbulent friction. 

![z_delta_max](img/z_delta_max.png)

The path to one of the neighbors (S_n) equals the path to the base (S_b) plus the path from base to the neighbor (S_bn) (Eq. 9). 

![S_n_eq1](img/S_n_eq1.png)

As there are many possibilities for the path from the starting point to our current base, we just take into account the shortest path, which corresponds to the highest Z_delta in the base. If Z_delta,max is set to infinity, or as in the code to 8898 m (= Mount Everest), we can calculate the shortest path from the starting point to our base with Eq. 10. 

![S_bn](img/S_bn.png)

With this equations we determine the maximum run out distance for the process, in the next steps we will explain how we handle and calculate the spreading on the 3D grid.



### Persistence Function

The persistence function P_i aims to reproduce the behavior of inertia, and weights the flow 
direction based on the change in direction with respect to the previous direction [3].
We introduced to scale the direction with Z_delta of the incoming direction (Z_delta,parent), 
so the direction from a cell with higher Z_delta will have more affect to the directions where the base cell spreads.

![](img/persistence.png)

The weights are defined by the cosine of the angle between parent, base and child/neighbor minus pi:  
![](img/persistence_cos_function.png)

So there are max. 3 childs that get input via the persistence function from one parent.

In the first calculation step, at the release or start cell, the persistence is set to one, because there exists no parent. So the first calculation step is totally dependent on the terrain.

### Terrain based routing

The terrain based routing is dependent on the slope angle. The Holmgren (1994) algorithm [1] is used in different kind of models and works well for avalanches but also rockfall or soil slides.
The exponent exp allows to control the divergence of the spreading. For avalanches a exponent of 8 shows good results.
To reach a single flow in step terrain (rockfall, soil slides, steepest descend), an exponent of 75 is considered.

![Holmgrem](img/flow_direction.png)
*Holmgrem Algorithm from 1994 [1]*



![Phi_Formula](img/Phi.png)



### Flux 

The values given by the terrain based routing and the persistence function are combined according to Eq.(16).

![](img/flux.png)

, where i is the direction and n are the neighbors from 1 to 8. R_i is then the flux that flows in direction i.
R_b is the flux in the base, for a release cell or starting cell the flux of the base equals one. \
The result of Eq. (16) is a 3 x 3 array with assigned flux values. A normalization stage is then 
required to bring the sum of the R_i's to the value of R_b. This aims at avoiding loss of flux [2].

### References

[1] [Holmgren, P. (1994).](https://www.researchgate.net/publication/229484151_Multiple_flow_direction_algorithms_for_runoff_modelling_in_grid_based_elevation_models_An_empirical_evaluation) 
Multiple flow direction algorithms for runoff modelling in
grid based elevation models: an empirical evaluation. Hydrological Processes, 8:327–334.


[2] [Horton, P., Jaboyedoff, M.,
Rudaz, B., and Zimmermann, M. (2013).](https://nhess.copernicus.org/articles/13/869/2013/nhess-13-869-2013.pdf) 
Flow-R, a model for susceptibility mapping of debris
flows and other gravitational hazards at a regional scale. Natural Hazards and Earth System
Science, 13:869–885.

[3] [Gamma, P. (1999).](https://www.researchgate.net/publication/34432465_dfwalk-Ein_Murgang-Simulationsprogramm_zur_Gefahrenzonierung) dfwalk - Ein
Murgang-Simulationsprogramm zur Gefahrenzonierung. PhD thesis, Universität Bern.

## Contact

For Questions contact:  
Michael Neuhauser, Austrian Research Centre for Forest: Michael.Neuhauser@bfw.gv.at  
Christopher D'Amboise, Austrian Research Centre for Forest: Christopher.DAmboise@bfw.gv.at