# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:44:52 2019

@author: Michael Neuhauser
This file is for reading and writing raster files


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
"""

import rasterio
import sys

def read_header(input_file):
    #Reads in the header of the raster file, input: filepath

    raster = rasterio.open(input_file)
    if raster is None:
        print('Unable to open {}'.format(input_file))
        sys.exit(1)

    header = {}
    header['ncols'] = raster.width
    header['nrows'] = raster.height
    header['xllcorner'] = (raster.transform * (0, 0))[0]
    header['yllcorner'] = (raster.transform * (0, raster.height))[1]
    header['cellsize'] = raster.transform[0]
    header['noDataValue'] = raster.nodata
    return header


def read_raster(input_file):

    header = read_header(input_file)
    raster = rasterio.open(input_file)
    my_array = raster.read(1)

    return my_array, header


def output_raster(file, file_out, raster):
    """Input is the original file, path to new file, raster_data
    
    Input parameters:
        file        the path to the file to reference on, mostly DEM on where 
                    Calculations were done
        file_out    path for the outputfile, possible extends are .asc or .tif"""

    raster_trans = rasterio.open(file)
    try:
        crs = rasterio.crs.CRS.from_dict(raster_trans.crs.data)
    except:
        crs = rasterio.crs.CRS.from_epsg(4326)
    if file_out[-3:] == 'asc': 
        new_dataset = rasterio.open(file_out, 'w', driver='AAIGrid', height = raster.shape[0], width = raster.shape[1], count=1,  dtype = raster.dtype, crs=crs, transform=raster_trans.transform, nodata=-9999)
    if file_out[-3:] == 'tif': 
        new_dataset = rasterio.open(file_out, 'w', driver='GTiff', height = raster.shape[0], width = raster.shape[1], count=1,  dtype = raster.dtype, crs=crs, transform=raster_trans.transform, nodata=-9999)
# =============================================================================
#     if file_out[-3:] != 'tif' & file_out[-3:] != 'asc':
#         print('This Fileformat is not supported: .{}'.format(file_out[-3:]))
#         return
# =============================================================================
    new_dataset.write(raster, 1)
    new_dataset.close()


#path = '/home/lawinenforschung/Desktop/PAR6_ValsGries_AUT/dhm_10_6.tif'
#output = 'example.asc'
#raster, header = read_raster(path)
#output_raster(path, output, raster)
#rasterio.rio.convert(output, 'example.asc', driver='AAIGrid')
