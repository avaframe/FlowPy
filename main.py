#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

@author: Neuhauser
"""
import sys 
import numpy as np
from datetime import datetime
from xml.etree import ElementTree as ET

import raster_io as io
import Simulation as Sim

from PyQt5.QtCore import pyqtSlot, QCoreApplication
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.QtGui import QIcon

FORM_CLASS = uic.loadUiType("Flow_GUI.ui")[0]


class GUI(QtWidgets.QMainWindow, FORM_CLASS):

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        # self.showMaximized()
        self.setWindowTitle("Flow Py GUI")
        self.setWindowIcon(QIcon('logo.jpg'))

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
        self.cpu_count = 1
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
        name = QFileDialog.getSaveFileName(self, 'Save File',
                                           "xml (*.xml);;All Files (*.*)")[0]

        root = ET.Element('root')
        wdir = ET.SubElement(root, 'wDir')
        dhm = ET.SubElement(root, 'DHM')
        release = ET.SubElement(root, 'Release')
        infra = ET.SubElement(root, 'Infrastrucutre')
        forest = ET.SubElement(root, 'Forest')

        wdir.set('Directory', 'Working')
        dhm.set('Directory', 'DHM')
        release.set('Directory', 'Release')
        infra.set('Directory', 'Infrastructure')
        forest.set('Directory', 'Forest')
        
        wdir.text = self.wDir_lineEdit.text()
        dhm.text = self.DEM_lineEdit.text()
        release.text = self.release_lineEdit.text()
        infra.text = self.infra_lineEdit.text()
        forest.text = self.forest_lineEdit.text()

        tree = ET.ElementTree(root)
        tree.write(name)

    def load(self):
        xml_file = QFileDialog.getOpenFileNames(self, 'Open xml',
                                                self.directory,
                                                "xml (*.xml);;All Files (*.*)")[0]
        
        tree = ET.parse(xml_file[0])
        root = tree.getroot()
        
        try:
            self.wDir_lineEdit.setText(root[0].text)
            self.directory = root[0].text
        except:
            print("No Working Directory Path in File!")
        
        try:
            self.DEM_lineEdit.setText(root[1].text)
        except:
            print("No DEM Path in File!")
            
        try:
            self.release_lineEdit.setText(root[2].text)
        except:
            print("No Release Path in File!")
            
        try:
            self.infra_lineEdit.setText(root[3].text)
        except:
            print("No Infrastructure Path in File!")
            
        try:
            self.forest_lineEdit.setText(root[4].text)
        except:
            print("No Forest Path in File!")

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
        if self.DEM_lineEdit.text() == '':
            self.showdialog('DEM Layer')
            return
        if self.release_lineEdit.text() == '':
            self.showdialog('Release Layer')
            return
        # Disable all input line Edits and Buttons
        self.calc_Button.setEnabled(False)
        self.wDir_lineEdit.setEnabled(False)
        self.DEM_lineEdit.setEnabled(False)
        self.release_lineEdit.setEnabled(False)
        self.infra_lineEdit.setEnabled(False)
        self.forest_lineEdit.setEnabled(False)
        

        # Start of Calculation
        # Read in raster files
        dem, header = io.read_raster(self.DEM_lineEdit.text())
        release, release_header = io.read_raster(self.release_lineEdit.text())
        
        #Check if Layers have same size!!!
        if header['ncols'] == release_header['ncols'] and header['nrows'] == release_header['nrows']:
            print("DEM and Release Layer ok!")
        else:
            print("Error: Release Layer doesn't match DEM!")
            return
                
        try:
            infra, infra_header = io.read_raster(self.infra_lineEdit.text())
            if (header['ncols'] == infra_header['ncols'] and header['nrows'] == infra_header['nrows']):
                print("DEM and Release ok!")
            else:
                print("Error: Infra Layer doesn't match DEM!")
                return
        except:
            infra = np.zeros_like(dem)
            
        try:
            forest, forest_header = io.read_raster(self.forest_lineEdit.text())
            if (header['ncols'] == forest_header['ncols'] and header['nrows'] == forest_header['nrows']):
                print("DEM and Release ok!")
            else:
                print("Error: Forest Layer doesn't match DEM!")
                return
        except:
            forest = np.zeros_like(dem)

        process = self.process_Box.currentText()
        self.elh = np.zeros_like(dem)
        self.mass = np.zeros_like(dem)
        self.cell_counts = np.zeros_like(dem)
        self.elh_sum = np.zeros_like(dem)

        # Calculation
        self.calc_class = Sim.Simulation(dem, header, release, release_header, forest, process)
        self.calc_class.value_changed.connect(self.update_progressBar)
        self.calc_class.finished.connect(self.thread_finished)
        self.calc_class.start()
                    
    def thread_finished(self, elh, mass, count_array, elh_sum):
        for i in range(len(elh)):
            self.elh = np.maximum(self.elh, elh[i])
            self.mass = np.maximum(self.mass, mass[i])
            self.cell_counts += count_array[i]
            self.elh_sum += elh_sum[i]
        self.output()
    
    def output(self):        
        if self.process_Box.currentText() == 'Avalanche':
            proc = 'ava'
        if self.process_Box.currentText() == 'Rockfall':
            proc = 'rf'
        if self.process_Box.currentText() == 'Soil Slides':
            proc = 'ds'
        time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/mass_{}_{}{}".format(proc, time_string, self.outputBox.currentText()), self.mass)
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/elh_{}_{}{}".format(proc, time_string, self.outputBox.currentText()), self.elh)
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/cell_counts_{}_{}{}".format(proc, time_string, self.outputBox.currentText()), self.cell_counts)
        io.output_raster(self.DEM_lineEdit.text(), self.directory + "/elh_sum_{}_{}{}".format(proc, time_string, self.outputBox.currentText()), self.elh_sum)
        print("Calculation finished")
        
        # Handle GUI
        self.progressBar.setValue(100)
        self.calc_Button.setEnabled(True)
        self.wDir_lineEdit.setEnabled(True)
        self.DEM_lineEdit.setEnabled(True)
        self.release_lineEdit.setEnabled(True)
        self.forest_lineEdit.setEnabled(True)


def main():
    """Go!"""
    app = QtWidgets.QApplication(sys.argv)
    ex = GUI()
    ex.show()
    app.exec_()


if __name__ == '__main__':
    main()
