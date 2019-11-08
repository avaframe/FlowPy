#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:41:45 2019

@author: Neuhauser
"""

import multiprocessing as mp
import gravi_core_gui as gc

from PyQt5.QtCore import QThread, pyqtSignal


class Simulation(QThread):
    value_changed = pyqtSignal(float)
    finished = pyqtSignal(list, list, list, list)

    def __init__(self, dem, header, release, release_header, forest, process):
        QThread.__init__(self)
        self.dem = dem
        self.header = header
        self.release = release
        self.release_header = release_header
        self.forest = forest
        self.process = process
        self.numberofprocesses = mp.cpu_count()

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
        
        # This part will is for Calculation of the top release cells and ereasing the lower ones
        if __name__ != '__main__':  # needed that it runs on windows
            
            release_list = gc.split_release(self.release, self.release_header)
            iterable = []
            for i in range(len(release_list)):
                iterable.append((self.dem, self.header, self.forest, self.process, release_list[i]))
        
            print("{} Processes started.".format(len(release_list)))
            pool = mp.Pool(len(release_list))
            results = pool.map(gc.calculation, iterable)
            pool.close()
            pool.join()
        
            print("Processes finished")
                
            elh_list = []
            mass_list = []
            cc_list = []
            elh_sum_list = []
            for i in range(len(results)):
                res = results[i]
                res = list(res)
                elh_list.append(res[0])
                mass_list.append(res[1])
                cc_list.append(res[2])
                elh_sum_list.append(res[3])
    
            self.finished.emit(elh_list, mass_list, cc_list, elh_sum_list)
            print("Results passed")
