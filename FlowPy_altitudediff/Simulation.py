#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:41:45 2019

@author: Michael Neuhauser
In this class the mutliprocessing is handled.


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
# Flow-Py libraries
import multiprocessing as mp
import flow_core as fc
import psutil
import sys

# PyQt libraries
from PyQt5.QtCore import QThread, pyqtSignal


class Simulation(QThread):
    value_changed = pyqtSignal(float)
    finished = pyqtSignal(list, list, list, list, list, list, list)

    def __init__(self, dem, header, release, release_header, infra, forest, calc_bool, alpha, exp, flux, max_z):
        QThread.__init__(self)
        self.dem = dem
        self.header = header
        self.release = release
        self.release_header = release_header
        self.infra = infra
        self.forest = forest
        self.numberofprocesses = mp.cpu_count()
        self.calc_bool = calc_bool
        self.alpha = alpha
        self.exp = exp
        self.flux = flux
        self.max_z = max_z

        avaiable_memory = psutil.virtual_memory()[1]
        dem_memory = sys.getsizeof(self.dem)
        release_memory = sys.getsizeof(self.release)
        infra_memory = sys.getsizeof(self.infra)
        puffer = avaiable_memory * 0.2  # ~ 5GB
        needed_memory = dem_memory * 7 + dem_memory + release_memory + infra_memory
        self.max_number_procces = int((avaiable_memory - puffer) / (needed_memory))
        print(avaiable_memory)

        print(
            "There are {} Bytes of Memory avaiable and {} Bytes needed per process. Max. Nr. of Processes = {}".format(
                avaiable_memory, needed_memory, self.max_number_procces))

    def run(self):




        # This part will is for Calculation of the top release cells and erasing the lower ones
        #if __name__ != '__main__':  # needed that it runs on windows, but it doesnt!!! if __name__ == main: would it be.
        if self.calc_bool:
            release_list = fc.split_release(self.release, self.release_header, min(mp.cpu_count() * 2, self.max_number_procces))  # mp.cpu_count() * 2,

            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(len(release_list))
            results = pool.map(fc.calculation, [[self.dem, self.header, self.infra, self.forest, release_pixel, self.alpha, self.exp, self.flux, self.max_z] for release_pixel in release_list])
            pool.close()
            pool.join()
        else:
            release_list = fc.split_release(self.release, self.release_header, min(mp.cpu_count() * 4, self.max_number_procces))  # mp.cpu_count() * 4,

            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(mp.cpu_count())
            #results = pool.map(gc.calculation, iterable)
            results = pool.map(fc.calculation_effect,
                               [[self.dem, self.header, self.forest, release_pixel, 
                                 self.alpha, self.exp, self.flux, self.max_z] 
                                for release_pixel in release_list])
            pool.close()
            pool.join()

        print("Processes finished")

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

        self.finished.emit(z_delta_list, flux_list, cc_list, z_delta_sum_list, backcalc_list, fp_ta_list, sl_ta_list)
        print("Results passed")
