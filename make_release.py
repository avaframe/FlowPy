# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:31:27 2021

@author: neuhauser
"""

import raster_io as io
import numpy as np

dem, header = io.read_raster("./examples/Oberammergau/PAR3_OAG_DGM_utm32n.tif")
slope, slope_header = io.read_raster("./examples/Oberammergau/slope.tif")

release = np.zeros_like(dem)

release[(dem > 1500) & (slope > 40) & (slope < 55)] = 1 

io.output_raster("./examples/Oberammergau/PAR3_OAG_DGM_utm32n.tif", "examples/Oberammergau/release.tif", release)