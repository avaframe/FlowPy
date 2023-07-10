#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

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
# import standard libraries
import os
import sys
import psutil
import shutil
import numpy as np
from datetime import datetime
from multiprocessing import cpu_count
import multiprocessing as mp
import logging
from xml.etree import ElementTree as ET
import pickle
import cProfile

# Flow-Py Libraries
import raster_io as io
#import Simulation as Sim
import flow_core as fc
import split_and_merge as SPAM


def main(args, kwargs):

    alpha = args[0]
    exp = args[1]
    directory = args[2]
    dem_path = args[3]
    release_path = args[4]
    if 'infra' in kwargs:
        infra_path = kwargs.get('infra')
        #print(infra_path)
    else:
        infra_path = None

    if 'flux' in kwargs:
        flux_threshold = float(kwargs.get('flux'))
        #print(flux_threshold)
    else:
        flux_threshold = 3 * 10 ** -4

    if 'max_z' in kwargs:
        max_z = kwargs.get('max_z')
        #print(max_z)
    else:
        max_z = 8848
        # Recomendet values:
                # Avalanche = 270
                # Rockfall = 50
                # Soil Slide = 12
    if 'maxCores' in kwargs:
        if kwargs.get('maxCores') == 'all':
            nCPU = mp.cpu_count()
        else:
            nCPU = mp.cpu_count()-1 #leaving 1 core "open" by default
    else:
        nCPU = mp.cpu_count()-1 #leaving 1 core "open" by default

    print("Starting...")
    print("...")

    start = datetime.now().replace(microsecond=0)
    infra_bool = False
    # Create result directory
    time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    try:
        os.makedirs(directory + 'res_{}/'.format(time_string))
        res_dir = ('res_{}/'.format(time_string))
    except FileExistsError:
        res_dir = ('res_{}/'.format(time_string))

    try:
        os.makedirs(directory + res_dir + 'temp/')
        temp_dir = (directory + res_dir + 'temp/')
    except FileExistsError:
        temp_dir = (directory + res_dir + 'temp/')

    # Setup logger
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=(directory + res_dir + 'log_{}.txt').format(time_string),
                        filemode='w')

    # Start of Calculation
    logging.info("--FlowPy development Version 'schuwa'--") #useful to utilized the Version of FlowPy
    logging.info('Start Calculation')
    logging.info('Alpha Angle: {}'.format(alpha))
    logging.info('Exponent: {}'.format(exp))
    logging.info('Flux Threshold: {}'.format(flux_threshold))
    logging.info('Max Z_delta: {}'.format(max_z))
    #logging.info
    # Read in raster files
    try:
        dem, header = io.read_raster(dem_path)
        logging.info('DEM File: {}'.format(dem_path))
    except FileNotFoundError:
        print("DEM: Wrong filepath or filename")
        return

    try:
        release, release_header = io.read_raster(release_path)
        logging.info('Release File: {}'.format(release_path))
    except FileNotFoundError:
        print("Wrong filepath or filename")
        return

    # Check if Layers have same size!!!
    if header['ncols'] == release_header['ncols'] and header['nrows'] == release_header['nrows']:
        print("DEM and Release Layer ok!")
    else:
        print("Error: Release Layer doesn't match DEM!")
        return

    try:
        infra, infra_header = io.read_raster(infra_path)
        if header['ncols'] == infra_header['ncols'] and header['nrows'] == infra_header['nrows']:
            print("Infra Layer ok!")
            infra_bool = True
            logging.info('Infrastructure File: {}'.format(infra_path))
        else:
            print("Error: Infra Layer doesn't match DEM!")
            return
    except:
        infra = np.zeros_like(dem)

    #del dem, infra
    logging.info('Files read in')

    cellsize = header["cellsize"]
    nodata = header["noDataValue"]
    tileCOLS = int(15000 / cellsize)
    tileROWS = int(15000 / cellsize)
    U = int(5000 / cellsize) # 5km overlap



    if (header['ncols'] * cellsize > 25000) or (header['nrows'] * cellsize > 25000):

        logging.info("Tiling is ON, because DEM dimensions larger than 15km in x and/or y")
        logging.info("Start Tiling.")
        print("Start Tiling...")

        SPAM.tileRaster(dem_path, "dem", temp_dir, tileCOLS, tileROWS, U)
        SPAM.tileRaster(release_path, "init", temp_dir, tileCOLS, tileROWS, U, isInit=True)
        if infra_bool:
            SPAM.tileRaster(infra_path, "infra", temp_dir, tileCOLS, tileROWS, U)

        print("Finished Tiling...")
        nTiles = pickle.load(open(temp_dir + "nTiles", "rb"))

        optList = []
        # das hier ist die batch-liste, die von mulitprocessing
        # abgearbeitet werden muss - sieht so aus:
        # [(0,0,alpha,exp,cellsize,-9999.),
        # (0,1,alpha,exp,cellsize,-9999.),
        # etc.]

        for i in range(nTiles[0]+1):
            for j in range(nTiles[1]+1):
                optList.append((i, j, alpha, exp, cellsize, nodata, flux_threshold, max_z, temp_dir))

        # Calculation
        logging.info('Multiprocessing starts, used cores: {}'.format(nCPU)) #this is only correct if number of tiles >= number of available cores!!
        print("{} Processes started and {} calculations to perform.".format(nCPU, len(optList)))
        pool = mp.Pool(nCPU)



        if infra_bool:
            pool.map(fc.calculation, optList)
        else:
            pool.map(fc.calculation_effect, optList)

        pool.close()
        pool.join()


        logging.info('Calculation finished, merging results.')

        # Merge calculated tiles
        z_delta = SPAM.MergeRaster(temp_dir, "res_z_delta")
        flux = SPAM.MergeRaster(temp_dir, "res_flux")
        cell_counts = SPAM.MergeRaster(temp_dir, "res_count")
        z_delta_sum = SPAM.MergeRaster(temp_dir, "res_z_delta_sum")
        fp_ta = SPAM.MergeRaster(temp_dir, "res_fp")
        sl_ta = SPAM.MergeRaster(temp_dir, "res_sl")
        if infra_bool:
            backcalc = SPAM.MergeRaster(temp_dir, "res_backcalc")

    else:
        print("No Tiling!")
        logging.info("No Tiling!")
        if infra_bool:
            release_list = fc.split_release(release, release_header, nCPU)

            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(len(release_list))
            results = pool.map(fc.calculation_small,
                               [[dem, header, infra, release_pixel, alpha, exp, flux_threshold, max_z]
                                for release_pixel in release_list])
            pool.close()
            pool.join()
        else:
            release_list = fc.split_release(release, release_header, nCPU)

            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(nCPU)
            # results = pool.map(gc.calculation, iterable)
            results = pool.map(fc.calculation_effect_small,
                               [[dem, header, release_pixel, alpha, exp, flux_threshold, max_z] for
                                release_pixel in release_list])
            pool.close()
            pool.join()

        z_delta = np.zeros_like(dem)
        flux = np.zeros_like(dem)
        cell_counts = np.zeros_like(dem)
        z_delta_sum = np.zeros_like(dem)
        backcalc = np.zeros_like(dem)
        fp_ta = np.zeros_like(dem)
        sl_ta = np.zeros_like(dem)
#
        z_delta_list = []
        flux_list = []
        cc_list = []
        z_delta_sum_list = []
        backcalc_list = []
        fp_ta_list = []
        sl_ta_list = []
        for i in range(len(results)):
            res = results[i]
            res = list(res)
            z_delta_list.append(res[0])
            flux_list.append(res[1])
            cc_list.append(res[2])
            z_delta_sum_list.append(res[3])
            backcalc_list.append(res[4])
            fp_ta_list.append(res[5])
            sl_ta_list.append(res[6])

        logging.info('Calculation finished, getting results.')
        for i in range(len(z_delta_list)):
            z_delta = np.maximum(z_delta, z_delta_list[i])
            flux = np.maximum(flux, flux_list[i])
            cell_counts += cc_list[i]
            z_delta_sum += z_delta_sum_list[i]
            backcalc = np.maximum(backcalc, backcalc_list[i])
            fp_ta = np.maximum(fp_ta, fp_ta_list[i])
            sl_ta = np.maximum(sl_ta, sl_ta_list[i])

    # time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    logging.info('Writing Output Files')
    output_format = '.tif'
    io.output_raster(dem_path,
                     directory + res_dir + "flux{}".format(output_format),
                     flux)
    io.output_raster(dem_path,
                     directory + res_dir + "z_delta{}".format(output_format),
                     z_delta)
    io.output_raster(dem_path,
                     directory + res_dir + "FP_travel_angle{}".format(output_format),
                     fp_ta)
    io.output_raster(dem_path,
                     directory + res_dir + "SL_travel_angle{}".format(output_format),
                     sl_ta)
    if not infra_bool:  # if no infra
        io.output_raster(dem_path,
                         directory + res_dir + "cell_counts{}".format(output_format),
                         cell_counts)
        io.output_raster(dem_path,
                         directory + res_dir + "z_delta_sum{}".format( output_format),
                         z_delta_sum)
    if infra_bool:  # if infra
        io.output_raster(dem_path,
                         directory + res_dir + "backcalculation{}".format(output_format),
                         backcalc)

    print("Calculation finished")
    print("...")
    end = datetime.now().replace(microsecond=0)
    logging.info('Calculation needed: ' + str(end - start) + ' seconds')
    logging.info('Deleting temp folder "%s" ...'%temp_dir)
    try:
        shutil.rmtree(temp_dir)
        logging.info('Deleted temp folder "%s"'%temp_dir)
    except OSError as e:
        print ("Error: %s : %s" %(temp_dir, e.strerror))


if __name__ == '__main__':
    #mp.set_start_method('spawn') # used in Windows
    argv = sys.argv[1:]
    #argv = ['--gui']
    #argv = ["25", "8", "./examples/dam/", "./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc", "./examples/dam/release_dam.tif"]
    #argv = ["15", "8", "./examples/dam/", "./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc", "./examples/dam/release_dam.tif", "infra=./examples/dam/infra.tif", "flux=0.0003", "max_z=270"]
    #argv = ["25", "8", "./examples/Arzler/", "./examples/Arzler/arzleralmdhm0101m_clipped.tif", "./examples/Arzler/release.tif"]
    #argv = ["25", "8", "./examples/Oberammergau/", "./examples/Oberammergau/PAR3_OAG_DGM_utm32n.tif", "./examples/Oberammergau/release.tif", "max_z=270"]
    #argv = ["25", "8", "./examples/Osttirol/", "./examples/Osttirol/DTM_5m.tif", "./examples/Osttirol/post_VAIA_release_areas_DGM_extend.tif", "max_z=270"]

    if len(argv) < 1:
    	print("Too few input arguments!!!")
    	sys.exit(1)
    if len(argv) == 1 and argv[0] == '--gui':
        Flow_Py_EXEC()
    else:
        args=[arg for arg in argv if arg.find('=')<0]
        kwargs={kw[0]:kw[1] for kw in [ar.split('=') for ar in argv if ar.find('=')>0]}

        #with cProfile.Profile() as cpr:
        main(args, kwargs)
        #cpr.dump_stats('~/profilingTest.pstats')

# example dam: python3 main.py 25 8 ./examples/dam/ ./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc ./examples/dam/release_dam.tif
# example dam: python3 main.py 25 8 ./examples/dam/ ./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc ./examples/dam/release_dam.tif infra=./examples/dam/infra.tif flux=0.0003 max_z=270
