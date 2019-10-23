#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

@author: Neuhauser
"""
import sys
import numpy as np

import raster_io as io
import gravi_core_gui as gc

from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox

FORM_CLASS = uic.loadUiType("Flow_GUI.ui")[0]


class GUI(QtWidgets.QMainWindow, FORM_CLASS):
    """extract the front and the Mask from a MTI Plot"""

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.showMaximized()
        self.setWindowTitle("Flow Py GUI")

        self.directory = '/home'

        self.alpha_Edit.setText('25')
        self.exp_Edit.setText('8')

        self.wDir_Button.clicked.connect(self.open_wDir)
        self.DEM_Button.clicked.connect(self.open_dhm)
        self.Release_Button.clicked.connect(self.open_release)
        self.infra_Button.clicked.connect(self.open_infra)
        self.forest_Button.clicked.connect(self.open_forest)
        self.process_Box.currentIndexChanged.connect(self.processChanged)
        self.calc_Button.clicked.connect(self.calculation)

    @pyqtSlot()
    def open_wDir(self):
        """Open the Working Directory, where results are stored"""
        self.directory = QFileDialog.getExistingDirectory(self, 'Open Working Directory',
                                                          '/home',
                                                          QFileDialog.ShowDirsOnly)
        self.wDir_lineEdit.setText(self.directory)

    def open_dhm(self):
        """Open the Working Directory, where results are stored"""
        dem_file = QFileDialog.getOpenFileNames(self, 'Open DEM',
                                                self.directory,
                                                "raster (*.asc);;All Files (*.*)")
        dem = dem_file[0]
        self.DEM_lineEdit.setText(dem[0])

    def open_release(self):
        """Open the Working Directory, where results are stored"""
        release_file = QFileDialog.getOpenFileNames(self, 'Open Release',
                                                    self.directory,
                                                    "raster (*.asc);;All Files (*.*)")
        release = release_file[0]
        self.release_lineEdit.setText(release[0])

    def open_infra(self):
        """Open the Working Directory, where results are stored"""
        infra_file = QFileDialog.getOpenFileNames(self, 'Open Infrastructure Layer',
                                                  self.directory,
                                                  "raster (*.asc);;All Files (*.*)")
        infra = infra_file[0]
        self.infra_lineEdit.setText(infra[0])

    def open_forest(self):
        """Open the Working Directory, where results are stored"""
        forest_file = QFileDialog.getOpenFileNames(self, 'Open Forest Layer',
                                                   self.directory,
                                                   "raster (*.asc);;All Files (*.*)")
        forest = forest_file[0]
        self.forest_lineEdit.setText(forest[0])

    def processChanged(self):
        if self.process_Box.currentText() == 'Avalanche':
            self.alpha_Edit.setText('25')
            self.exp_Edit.setText('8')
        if self.process_Box.currentText() == 'Rockfall':
            self.alpha_Edit.setText('32')
            self.exp_Edit.setText('75')
        if self.process_Box.currentText() == 'Soil Slides':
            self.alpha_Edit.setText('22')
            self.exp_Edit.setText('75')

    def showdialog(self, path):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("No " + path + " set")
        msg.setWindowTitle("Error")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def calculation(self):
        # Check if input is ok
        if self.wDir_lineEdit.text() == '':
            self.showdialog('Working Directory')
            return
        else:
            path = self.wDir_lineEdit.text()
        if self.DEM_lineEdit.text() == '':
            self.showdialog('DEM Layer')
            return
        else:
            dem_file = self.DEM_lineEdit.text()
        if self.release_lineEdit.text() == '':
            self.showdialog('Release Layer')
            return
        else:
            release_file = self.release_lineEdit.text()

        # Start of Calculation
        # Read in raster files
        dem, header = io.read_raster(self.DEM_lineEdit.text())
        cellsize = header["cellsize"]
        nodata = header["noDataValue"]
        release, header_release = io.read_raster(self.release_lineEdit.text())

        # infra, header = io.read_raster(infra_path) needed for backcalculation
        try:
            forest, header_forest = io.read_raster(self.forest_lineEdit.text())
        except:
            forest = np.zeros_like(dem)

        # Calculation
        # Create needed arrays
        mass_array = np.zeros_like(dem)
        elh = np.zeros_like(dem)
        count_array = np.zeros_like(dem)

        row_list, col_list = gc.get_start_idx(dem, release)
        startcell_idx = 0
        while startcell_idx < len(row_list):
            self.progressBar.value(round((startcell_idx + 1) / len(row_list) * 100, 2))
            elh, mass_array, count_array = gc.calculation(dem, release, forest, self.process_Box.currentText(), cellsize, nodata, row_list, col_list, startcell_idx, elh, mass_array, count_array)
            release[elh > 0] = 0  # Check if i hited a release Cell, if so set it to zero and get again the indexes of release cells
            # ToDo: if i hit a startcell add this "mass"
            # ToDo: Backcalulation
            row_list, col_list = gc.get_start_idx(dem, release)
            startcell_idx += 1

        # Output
        io.output_raster(dem_file, self.directory + "mass.tif", mass_array)
        io.output_raster(dem_file, self.directory + "elh.tif", elh)
        io.output_raster(dem_file, self.directory + "cell_count.tif", count_array)


def main():
    """Go!"""
    app = QtWidgets.QApplication(sys.argv)
    ex = GUI()
    ex.show()
    app.exec_()


if __name__ == '__main__':
    main()
