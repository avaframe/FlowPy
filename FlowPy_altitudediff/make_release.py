# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:31:27 2021

@author: neuhauser
"""

import raster_io as io
import numpy as np

dem, header = io.read_raster("./examples/ValsGries/PAR6_Vals_Gries_dtm_10_utm32n_bil_.tif")
slope, slope_header = io.read_raster("./examples/ValsGries/slope.asc")

release = np.zeros_like(dem)

release[(dem > 1800) & (slope > 30) & (slope < 35)] = 1 

io.output_raster("./examples/ValsGries/PAR6_Vals_Gries_dtm_10_utm32n_bil_.tif", "release.tif", release)