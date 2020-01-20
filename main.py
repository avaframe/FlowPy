#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

@author: Neuhauser
"""
# import standard libraries
import os
import sys
import numpy as np
from datetime import datetime
from multiprocessing import cpu_count
import logging
from xml.etree import ElementTree as ET

# Flow-Py Libraries
import raster_io as io
import Simulation as Sim

# Libraries for GUI, PyQt5
from PyQt5 import uic
from PyQt5.QtCore import pyqtSlot, QCoreApplication
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QMainWindow, QApplication
from PyQt5.QtGui import QIcon

FORM_CLASS = uic.loadUiType("Flow_GUI.ui")[0]


class GUI(QMainWindow, FORM_CLASS):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
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
        self.prot_for_bool = False
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
            
    def set_gui_bool(self, bool):
        self.calc_Button.setEnabled(bool)
        self.wDir_lineEdit.setEnabled(bool)
        self.DEM_lineEdit.setEnabled(bool)
        self.release_lineEdit.setEnabled(bool)
        self.infra_lineEdit.setEnabled(bool)
        self.forest_lineEdit.setEnabled(bool)

    @pyqtSlot()
    def save(self):
        """Save the input paths"""
        name = QFileDialog.getSaveFileName(self, 'Save File',
                                           ".xml")[0]

        root = ET.Element('root')
        wdir = ET.SubElement(root, 'wDir')
        dhm = ET.SubElement(root, 'DHM')
        release = ET.SubElement(root, 'Release')
        infra = ET.SubElement(root, 'Infrastructure')
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
        """Open digital elevation model"""
        dem_file = QFileDialog.getOpenFileNames(self, 'Open DEM',
                                                self.directory,
                                                "tif (*.tif);;raster (*.asc);;All Files (*.*)")
        dem = dem_file[0]
        self.DEM_lineEdit.setText(dem[0])

    def open_release(self):
        """Open release layer"""
        release_file = QFileDialog.getOpenFileNames(self, 'Open Release',
                                                    self.directory,
                                                    "tif (*.tif);;raster (*.asc);;All Files (*.*)")
        release = release_file[0]
        self.release_lineEdit.setText(release[0])

    def open_infra(self):
        """Open infrastructure layer"""
        infra_file = QFileDialog.getOpenFileNames(self, 'Open Infrastructure Layer',
                                                  self.directory,
                                                  "tif (*.tif);;raster (*.asc);;All Files (*.*)")
        infra = infra_file[0]
        self.infra_lineEdit.setText(infra[0])

    def open_forest(self):
        """Open forest layer"""
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

        self.progress_value = sum(self.thread_list) / len(self.thread_list)
        self.progressBar.setValue(self.progress_value)
        for i in range(len(self.thread_list)):
            sys.stdout.write(
                "Thread {}: Startcell {} of {} = {}%"'\n'.format(i + 1, self.start_list[i], self.end_list[i],
                                                                 self.thread_list[i]))
            sys.stdout.flush()
        for i in range(len(self.thread_list)):
            sys.stdout.write('\x1b[1A' + '\x1b[2K')  # Go 1 line up and erase that line

    @staticmethod
    def showdialog(path):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("No " + path + " set")
        msg.setWindowTitle("Error")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def calculation(self):
        self.start = datetime.now().replace(microsecond=0)
        calc_bool = False

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
        self.set_gui_bool(False)
# =============================================================================
#         self.calc_Button.setEnabled(False)
#         self.wDir_lineEdit.setEnabled(False)
#         self.DEM_lineEdit.setEnabled(False)
#         self.release_lineEdit.setEnabled(False)
#         self.infra_lineEdit.setEnabled(False)
#         self.forest_lineEdit.setEnabled(False)
# =============================================================================
        
                # Create result directory
        time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        try:
            os.makedirs(self.wDir_lineEdit.text() + '/res_{}/'.format(time_string))
            self.res_dir = ('/res_{}/'.format(time_string))
        except FileExistsError:
            self.res_dir = ('/res_{}/'.format(time_string))

            # Setup logger

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            filename=(self.directory + self.res_dir + 'log_{}.txt').format(time_string),
                            filemode='w')

        # Start of Calculation
        logging.info('Start Calculation')
        # Read in raster files
        try:
            dem, header = io.read_raster(self.DEM_lineEdit.text())
            logging.info('DEM File: {}'.format(self.DEM_lineEdit.text()))
        except FileNotFoundError:
            print("Wrong filepath or filename")
            self.set_gui_bool(True)
            return

        try:
            release, release_header = io.read_raster(self.release_lineEdit.text())
            logging.info('Release File: {}'.format(self.release_lineEdit.text()))
        except FileNotFoundError:
            print("Wrong filepath or filename")
            self.set_gui_bool(True)
            return

        # Check if Layers have same size!!!
        if header['ncols'] == release_header['ncols'] and header['nrows'] == release_header['nrows']:
            print("DEM and Release Layer ok!")
        else:
            print("Error: Release Layer doesn't match DEM!")
            self.set_gui_bool(True)
            return

        try:
            infra, infra_header = io.read_raster(self.infra_lineEdit.text())
            if header['ncols'] == infra_header['ncols'] and header['nrows'] == infra_header['nrows']:
                print("Infra Layer ok!")
                calc_bool = True
                logging.info('Infrastructure File: {}'.format(self.infra_lineEdit.text()))
            else:
                print("Error: Infra Layer doesn't match DEM!")
                self.set_gui_bool(True)
                return
        except:
            infra = np.zeros_like(dem)

        try:
            forest, forest_header = io.read_raster(self.forest_lineEdit.text())
            if header['ncols'] == forest_header['ncols'] and header['nrows'] == forest_header['nrows']:
                print("Forest Layer ok!")
                self.prot_for_bool = True
                logging.info('Forest File: {}'.format(self.forest_lineEdit.text()))
            else:
                print("Error: Forest Layer doesn't match DEM!")
                self.set_gui_bool(True)
                return
        except:
            forest = np.zeros_like(dem)


        logging.info('Files read in')

        process = self.process_Box.currentText()
        logging.info('Process: {}'.format(process))
        self.elh = np.zeros_like(dem)
        self.mass = np.zeros_like(dem)
        self.cell_counts = np.zeros_like(dem)
        self.elh_sum = np.zeros_like(dem)
        self.backcalc = np.zeros_like(dem)

        # Calculation
        self.calc_class = Sim.Simulation(dem, header, release, release_header, infra, forest, process, calc_bool)
        self.calc_class.value_changed.connect(self.update_progressBar)
        self.calc_class.finished.connect(self.thread_finished)
        logging.info('Multiprocessing starts, used cores: {}'.format(cpu_count()))
        self.calc_class.start()

    def thread_finished(self, elh, mass, count_array, elh_sum, backcalc):
        logging.info('Calculation finished, getting results.')
        for i in range(len(elh)):
            self.elh = np.maximum(self.elh, elh[i])
            self.mass = np.maximum(self.mass, mass[i])
            self.cell_counts += count_array[i]
            self.elh_sum += elh_sum[i]
            self.backcalc = np.maximum(self.backcalc, backcalc[i])
        self.output()

    def output(self):
        if self.process_Box.currentText() == 'Avalanche':
            proc = 'ava'
        if self.process_Box.currentText() == 'Rockfall':
            proc = 'rf'
        if self.process_Box.currentText() == 'Soil Slides':
            proc = 'ds'
        # time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        logging.info('Writing Output Files')
        io.output_raster(self.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "mass_{}{}".format(proc, self.outputBox.currentText()),
                         self.mass)
        io.output_raster(self.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "elh_{}{}".format(proc, self.outputBox.currentText()),
                         self.elh)
        io.output_raster(self.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "cell_counts_{}{}".format(proc, self.outputBox.currentText()),
                         self.cell_counts)
        io.output_raster(self.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "elh_sum_{}{}".format(proc, self.outputBox.currentText()),
                         self.elh_sum)
        io.output_raster(self.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "backcalculation_{}{}".format(proc, self.outputBox.currentText()),
                         self.backcalc)

        # Output of Protection forest if Infra structure and forest layer are
        # provided

        print("Calculation finished")
        end = datetime.now().replace(microsecond=0)
        logging.info('Calculation needed: ' + str(end - self.start) + ' seconds')

        # Handle GUI
        self.progressBar.setValue(100)
        self.set_gui_bool(True)
# =============================================================================
#         self.calc_Button.setEnabled(True)
#         self.wDir_lineEdit.setEnabled(True)
#         self.DEM_lineEdit.setEnabled(True)
#         self.release_lineEdit.setEnabled(True)
#         self.infra_lineEdit.setEnabled(True)
#         self.forest_lineEdit.setEnabled(True)
# =============================================================================


def main():
    """Go!"""
    app = QApplication(sys.argv)
    ex = GUI()
    ex.show()
    app.exec_()


if __name__ == '__main__':
    main()
