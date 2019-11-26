#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:41:45 2019

@author: Neuhauser
"""
# Flow-Py libraries
import multiprocessing as mp
import gravi_core_gui as gc

# PyQt libraries
from PyQt5.QtCore import QThread, pyqtSignal


class Simulation(QThread):
    value_changed = pyqtSignal(float)
    finished = pyqtSignal(list, list, list, list, list)

    def __init__(self, dem, header, release, release_header, infra, forest, process, calc_bool):
        QThread.__init__(self)
        self.dem = dem
        self.header = header
        self.release = release
        self.release_header = release_header
        self.infra = infra
        self.forest = forest
        self.process = process
        self.numberofprocesses = mp.cpu_count()
        self.calc_bool = calc_bool

    def run(self):

        # This part is for Calculation of all release cells
# =============================================================================
#         row_list, col_list = get_start_idx(self.dem, self.release)
#         divided_rowlist = list(divide_chunks(row_list, int(len(row_list)/self.numberofprocesses - 1)))
#         divided_collist = list(divide_chunks(col_list, int(len(col_list)/self.numberofprocesses - 1)))
# 
#         iterable = []
#         for i in range(self.numberofprocesses):
#             iterable.append((self.dem, self.header, self.forest, self.process, divided_rowlist[i], divided_collist[i]))
# =============================================================================
        
        # This part will is for Calculation of the top release cells and erasing the lower ones
        #if __name__ != '__main__':  # needed that it runs on windows, but it doesnt!!! if __name__ == main: would it be.
        if self.calc_bool:    
            release_list = gc.split_release(self.release, self.release_header, mp.cpu_count()*2)
    
            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(len(release_list))
            results = pool.map(gc.calculation, [[self.dem, self.header, self.infra, self.forest, self.process, release_pixel] for release_pixel in release_list])
            pool.close()
            pool.join()
        else:
            release_list = gc.split_release(self.release, self.release_header, mp.cpu_count()*4)
    
            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(mp.cpu_count())
            #results = pool.map(gc.calculation, iterable)
            results = pool.map(gc.calculation_effect, [[self.dem, self.header, self.forest, self.process, release_pixel] for release_pixel in release_list])
            pool.close()
            pool.join()

        print("Processes finished")

        elh_list = []
        mass_list = []
        cc_list = []
        elh_sum_list = []
        backcalc_list = []
        for i in range(len(results)):
            res = results[i]
            res = list(res)
            elh_list.append(res[0])
            mass_list.append(res[1])
            cc_list.append(res[2])
            elh_sum_list.append(res[3])
            backcalc_list.append(res[4])

        self.finished.emit(elh_list, mass_list, cc_list, elh_sum_list, backcalc_list)
        print("Results passed")
