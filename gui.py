#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

@author: Neuhauser
"""
import sys
import numpy as np
import multiprocessing
from multiprocessing import Pool
from xml.etree import ElementTree as ET

import raster_io as io
import gravi_core_gui as gc

from PyQt5.QtCore import pyqtSlot, QCoreApplication, QThreadPool
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox

FORM_CLASS = uic.loadUiType("Flow_GUI.ui")[0]


class GUI(QtWidgets.QMainWindow, FORM_CLASS):
    """extract the front and the Mask from a MTI Plot"""

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        #self.showMaximized()
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
        self.actionSave.triggered.connect(self.save)
        self.actionLoad.triggered.connect(self.load)
        self.actionQuit.triggered.connect(self.quit)

        self.calc_class = None
        self.threads_calc = 0
        self.progress_value = 0
        self.cpu_count = multiprocessing.cpu_count()
        self.thread_list = []
        self.start_list = []
        self.end_list = []
        for i in range(self.cpu_count):
            self.thread_list.append(0)
            self.start_list.append(0)
            self.end_list.append(0)

    @pyqtSlot()
    def save(self):
        """Save the input paths"""
        name = QFileDialog.getSaveFileName(self, 'Save File')[0]
        #file = open(name, 'w')
        #text = self.wDirlineEdit.text()

        root = ET.Element("root")
        wdir = ET.SubElement(root, "wDir")
        dhm = ET.SubElement(root, "DHM")

        ET.SubElement(wdir, self.wDir_lineEdit.text())
        ET.SubElement(dhm, self.DEM_lineEdit.text())

        tree = ET.ElementTree(root)
        tree.write(name)
        #text = self.textEdit.toPlainText()
        #file.write(text)
        #file.close()

    def load(self):
        xml_file = QFileDialog.getOpenFileNames(self, 'Open xml',
                                                self.directory,
                                                "xml (*.xml);;All Files (*.*)")[0]
        xml = xml_file[0]
        print(xml)

    def quit(self):
        QCoreApplication.quit()

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
                                                "tif (*.tif);;raster (*.asc);;All Files (*.*)")
        dem = dem_file[0]
        self.DEM_lineEdit.setText(dem[0])

    def open_release(self):
        """Open the Working Directory, where results are stored"""
        release_file = QFileDialog.getOpenFileNames(self, 'Open Release',
                                                    self.directory,
                                                    "tif (*.tif);;raster (*.asc);;All Files (*.*)")
        release = release_file[0]
        self.release_lineEdit.setText(release[0])

    def open_infra(self):
        """Open the Working Directory, where results are stored"""
        infra_file = QFileDialog.getOpenFileNames(self, 'Open Infrastructure Layer',
                                                  self.directory,
                                                  "tif (*.tif);;raster (*.asc);;All Files (*.*)")
        infra = infra_file[0]
        self.infra_lineEdit.setText(infra[0])

    def open_forest(self):
        """Open the Working Directory, where results are stored"""
        forest_file = QFileDialog.getOpenFileNames(self, 'Open Forest Layer',
                                                   self.directory,
                                                   "tif (*.tif);;raster (*.asc);;All Files (*.*)")
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

    def update_progressBar(self, float, thread, start, end):
        #print(float)
        self.thread_list[thread] = float
        self.start_list[thread] = start
        self.end_list[thread] = end
                    
        self.progress_value = sum(self.thread_list)/len(self.thread_list)
        self.progressBar.setValue(self.progress_value)
        for i in range(len(self.thread_list)):             
                sys.stdout.write("Thread {}: Startcell {} of {} = {}%"'\n'.format(i+1, self.start_list[i], self.end_list[i], self.thread_list[i]))
                sys.stdout.flush()
        for i in range(len(self.thread_list)):
            sys.stdout.write('\x1b[1A' + '\x1b[2K') # Go 1 line up and erase that line
                
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
        release, header_release = io.read_raster(self.release_lineEdit.text())
        # infra, header = io.read_raster(infra_path) needed for backcalculation
        try:
            forest, header_forest = io.read_raster(self.forest_lineEdit.text())
        except:
            forest = np.zeros_like(dem)
        process = self.process_Box.currentText()
        self.elh = np.zeros_like(dem)
        self.mass = np.zeros_like(dem)
        self.cell_counts = np.zeros_like(dem)

        # Calculation
        
        #for i in range(2):
# =============================================================================
#         self.calc_class = gc.Simulation(dem, header, release, forest, process)
#         self.calc_class.value_changed.connect(self.update_progressBar)
#         self.calc_class.finished.connect(self.output)
#         self.calc_class.start()
# =============================================================================
        
        # Try to set up multiple threads
        #cpu_count = multiprocessing.cpu_count()
        #pool = Pool(processes=cpu_count)
        self.threadpool = QThreadPool()
        thread_list = []
        for i in range(self.cpu_count):
            #proc = multiprocessing.Process(target = gc.Simulation, args=(dem, header, release, forest, process, i, cpu_count))
            thread_list.append(gc.Simulation(dem, header, release, forest, process, i, self.cpu_count))
            #thread_list.append(proc)
        for thread in thread_list:
            thread.signals.value_changed.connect(self.update_progressBar)
            thread.signals.finished.connect(self.thread_finished)
            self.threadpool.start(thread)
            #pool.map(thread)

        #self.threadpool.waitForDone(self.output)
        # Output
        #self.calc_class.finished.connect(self.output)
    
    def thread_finished(self, elh, mass_array, count_array):
        self.threads_calc +=1
        self.elh = np.maximum(self.elh, elh)
        self.mass = np.maximum(self.mass, mass_array)
        self.cell_counts += count_array 
        if self.threads_calc == self.cpu_count:
            self.output()
    
    def output(self):
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/mass_gui.tif", self.mass)
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/elh_gui.tif", self.elh)
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/cell_count_gui.tif", self.cell_counts)
        self.progressBar.setValue(100)
        print("Calculation finished")


def main():
    """Go!"""
    app = QtWidgets.QApplication(sys.argv)
    ex = GUI()
    ex.show()
    app.exec_()


if __name__ == '__main__':
    main()
