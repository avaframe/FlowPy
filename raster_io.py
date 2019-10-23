# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:44:52 2019

@author: Michael Neuhauser
This file is for reading and writing raster files
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
    """Input is the original file, path to new file, raster_data, and the EPSG Code"""

    raster_trans = rasterio.open(file)
    crs = rasterio.crs.CRS.from_dict(raster_trans.crs.data)
    new_dataset = rasterio.open(file_out, 'w', driver='GTiff', height = raster.shape[0], width = raster.shape[1], count=1,  dtype = raster.dtype, crs=crs, transform=raster_trans.transform, nodata=-9999)
    new_dataset.write(raster, 1)
    new_dataset.close()  

    
def output_raster_v2(file_out, raster, epsg):
    """Input is path to file, raster_data"""
    

    #raster_trans = rasterio.open(file)
    crs = rasterio.crs.CRS.from_epsg(epsg)
    new_dataset = rasterio.open(file_out, 'w', driver='GTiff', height = raster.shape[0], width = raster.shape[1], count=1,  dtype = raster.dtype, crs=crs, nodata=-9999)
    new_dataset.write(raster, 1)
    new_dataset.close()  


#header = read_header('/home/P/Projekte/18130-GreenRisk4Alps/Simulation/PAR3_Oberammergau/DEM_1_3.asc')
#output_raster_v2('example.tif', raster)