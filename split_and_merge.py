#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import logging
import pickle
import gc
import numpy as np
from dakumo_base_modules import IO_functionality as iof

# import tkFileDialog
# wDir = tkFileDialog.askdirectory(title="define working directory")


def tileRaster(fNameIn, fNameOut, dirName, xDim, yDim, U, isInit=False):

    if not os.path.exists(dirName):
        os.makedirs(dirName)

    largeRaster, largeHeader = iof.read_raster(fNameIn)
    # einlesen des Rasters und der Header


    i, j, imax, jmax = 0, 0, 0, 0
    sX, sY, eX, eY = 0, 0, 0, 0
    # starte mit den tiles in der NW-Ecke sX,eX = cols; sY,eY = rows;
    # cs=largeHeader['cellsize']
    # xllc = largeHeader['xllcorner']
    # yllc = largeHeader['yllcorner']

    if not largeHeader['noDataValue']:
        print('Error: NoDataValue for file {} is None!'.format(fNameIn))
        # Check for NoDataValue, raises error when merging

    nrows, ncols = largeRaster.shape[0], largeRaster.shape[1]
    if nrows != largeHeader['nrows']:
        print("Size of header rows and file rows do not match. File: {}".format(fNameIn))
        print("Header: {}, File: {}".format(largeHeader['nrows'], nrows))
    elif ncols != largeHeader['ncols']:
        print("Size of header columns and file columns do not match. File: {}".format(fNameIn))
        print("Header: {}, File: {}".format(largeHeader['ncols'], ncols))

    pickle.dump((nrows, ncols), open("%s/extentLarge" % (dirName), "wb"))

    # print ncols, nrows

    I, J, IMAX, JMAX = 0, 0, 0, 0

    while eY < nrows:
        eY = sY + yDim
        while eX < ncols:
            eX = sX + xDim

    # rangeRowsCols = ((sY,eY),(sX,eX))
    # pickle.dump(rangeRowsCols, open("%s/ext_%i_%i"%(dirName,i,j),"wb"))

    # headerTile = {}
    # headerTile['ncols'] = eX-sX
    # headerTile['nrows'] = eY-sY
    # headerTile['xllcorner'] = xllc + sX*cs
    # headerTile['yllcorner'] = yllc + nrows*cs - eY*cs
    # headerTile['cellsize'] = cs
    # headerTile['noDataValue'] = largeHeader['noDataValue']

    # pickle.dump( headerTile, open( "temp/header%d_%d.p"%(i,j), "wb" ) )
    # np.save("%s/%s_%i_%i"%(dirName,fNameOut, i, j), largeRaster[sY:eY,sX:eX])
    # pickle.dump(, open( "temp/header_large.p"%(fNameOut, i,j), "wb" ) )
    # logging.info("saved %s - TileNr.: %i_%i"%(fNameOut,i,j))

            sX = eX-2*U
            JMAX = max(J, JMAX)
            J += 1
        sX, J, eX = 0, 0, 0
        sY = eY-2*U
        IMAX = max(I, IMAX)
        I += 1

    sX, sY, eX, eY = 0, 0, 0, 0

    if isInit is False:
        while eY < nrows:
            eY = sY+yDim
            while eX < ncols:
                eX = sX+xDim
                rangeRowsCols = ((sY, eY), (sX, eX))
                pickle.dump(rangeRowsCols,
                            open("%s/ext_%i_%i" % (dirName, i, j), "wb"))

                # headerTile = {}
                # headerTile['ncols'] = eX-sX
                # headerTile['nrows'] = eY-sY
                # headerTile['xllcorner'] = xllc + sX*cs
                # headerTile['yllcorner'] = yllc + nrows*cs - eY*cs
                # headerTile['cellsize'] = cs
                # headerTile['noDataValue'] = largeHeader['noDataValue']

                # pickle.dump( headerTile,
                # open( "temp/header%d_%d.p"%(i,j), "wb" ) )
                np.save("%s/%s_%i_%i" % (dirName, fNameOut, i, j),
                        largeRaster[sY:eY, sX:eX])
                # pickle.dump(,
                # open( "temp/header_large.p"%(fNameOut, i,j), "wb" ) )
                logging.info("saved %s - TileNr.: %i_%i", fNameOut, i, j)

                sX = eX-2*U
                jmax = max(j, jmax)
                j += 1
            sX, j, eX = 0, 0, 0
            sY = eY-2*U
            imax = max(i, imax)
            i += 1
    else:
        while eY < nrows:
            eY = sY+yDim
            while eX < ncols:
                eX = sX+xDim

                rangeRowsCols = ((sY, eY), (sX, eX))
                pickle.dump(rangeRowsCols,
                            open("%s/ext_%i_%i" % (dirName, i, j), "wb"))

                # headerTile = {}
                # headerTile['ncols'] = eX-sX
                # headerTile['nrows'] = eY-sY
                # headerTile['xllcorner'] = xllc + sX*cs
                # headerTile['yllcorner'] = yllc + nrows*cs - eY*cs
                # headerTile['cellsize'] = cs
                # headerTile['noDataValue'] = largeHeader['noDataValue']
                # pickle.dump(headerTile,
                # open( "temp/hd_%s%_d_%d.p"%(fNameOut, i,j), "wb" ) )

                initRas = largeRaster[sY:eY, sX:eX].copy()
                # shapeX = np.shape(initRas)[1]
                # shapeY = np.shape(initRas)[0]
                if j != JMAX:
                    initRas[:, -U:] = -9999  # Rand im Osten
                if i != 0:
                    initRas[0:U, :] = -9999  # Rand im Norden
                if j != 0:
                    initRas[:, 0:U] = -9999  # Rand im Westen
                if i != IMAX:
                    initRas[-U:, :] = -9999  # Rand im Sueden

                # logging.info("%i_%i"%(shapeX-U, shapeX))
                np.save("%s/%s_%i_%i" % (dirName, fNameOut, i, j), initRas)
                del initRas
                # pickle.dump(,
                # open( "temp/header_large.p"%(fNameOut, i,j), "wb" ) )
                logging.info("saved %s - TileNr.: %i_%i", fNameOut, i, j)

                sX = eX-2*U
                jmax = max(j, jmax)
                j += 1
            sX, j, eX = 0, 0, 0
            sY = eY-2*U
            imax = max(i, imax)
            i += 1

    pickle.dump((imax, jmax), open("%s/nTiles" % (dirName), "wb"))
    logging.info("finished tiling %s: nTiles=%s\n----------------------------",
                 fNameIn, (imax+1)*(jmax+1))

    # del largeRaster, largeHeader
    del largeHeader, largeRaster
    gc.collect()
    # return largeRaster


def MergeRaster(inDirPath, fName):

    os.chdir(inDirPath)

    extL = pickle.load(open("extentLarge", "rb"))
    # print extL
    nTiles = pickle.load(open("nTiles", "rb"))

    mergedRas = np.zeros((extL[0], extL[1]))
    # create Raster with original size
    mergedRas[:, :] = np.NaN

    for i in range(nTiles[0]+1):
        for j in range(nTiles[1]+1):
            smallRas = np.load("%s_%i_%i.npy" % (fName, i, j))
            # print smallRas
            pos = pickle.load(open("ext_%i_%i" % (i, j), "rb"))
            # print pos

            mergedRas[pos[0][0]:pos[0][1], pos[1][0]:pos[1][1]] =\
                np.fmax(mergedRas[pos[0][0]:pos[0][1],
                                  pos[1][0]:pos[1][1]], smallRas)
            del smallRas
            logging.info("appended result %s_%i_%i", fName, i, j)

    return mergedRas
    del mergedRas
