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

# PyQt libraries
from PyQt5.QtCore import QThread, pyqtSignal


class Simulation(QThread):
    value_changed = pyqtSignal(float)
    finished = pyqtSignal()

    def __init__(self, optList, infra_bool):
        QThread.__init__(self)
        self.optList = optList
        self.infra_bool = infra_bool


    def run(self):

        print("{} Processes started.".format(mp.cpu_count() - 1))
        pool = mp.Pool(mp.cpu_count() - 1)
        pool.map(fc.calculation, self.optList)
        pool.close()
        pool.join()

        print("Processes finished")
        self.finished.emit()