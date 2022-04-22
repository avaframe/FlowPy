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
import numpy as np
from datetime import datetime
from multiprocessing import cpu_count
import multiprocessing as mp
import logging
from xml.etree import ElementTree as ET

# Flow-Py Libraries
import raster_io as io
import Simulation as Sim
import flow_core as fc

# Libraries for GUI, PyQt5
from PyQt5.QtCore import pyqtSlot, QCoreApplication
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QMainWindow, QApplication

from Flow_GUI import Ui_MainWindow


class Flow_Py_EXEC():

    def __init__(self):
        
        app = QApplication(sys.argv) 
        MainWindow = QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(MainWindow)
        
        # self.showMaximized()
        #self.ui.setWindowTitle("Flow-Py")
        #self.ui.setWindowIcon(QIcon('logo.jpg'))

        self.directory = os.getcwd()

        self.ui.alpha_Edit.setText('25')
        self.ui.exp_Edit.setText('8')
        self.ui.flux_Edit.setText('0.003')
        self.ui.z_Edit.setText('8848')

        self.ui.wDir_Button.clicked.connect(self.open_wDir)
        self.ui.DEM_Button.clicked.connect(self.open_dhm)
        self.ui.Release_Button.clicked.connect(self.open_release)
        self.ui.infra_Button.clicked.connect(self.open_infra)
        self.ui.forest_Button.clicked.connect(self.open_forest)
        #self.ui.process_Box.currentIndexChanged.connect(self.processChanged)
        self.ui.calc_Button.clicked.connect(self.calculation)
        self.ui.actionSave.triggered.connect(self.save)
        self.ui.actionLoad.triggered.connect(self.load)
        self.ui.actionQuit.triggered.connect(self.quit)

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
            
        # show the constructed window
        MainWindow.show()
        sys.exit(app.exec_())
            
    def set_gui_bool(self, bool):
        self.ui.calc_Button.setEnabled(bool)
        self.ui.wDir_lineEdit.setEnabled(bool)
        self.ui.DEM_lineEdit.setEnabled(bool)
        self.ui.release_lineEdit.setEnabled(bool)
        self.ui.infra_lineEdit.setEnabled(bool)
        self.ui.forest_lineEdit.setEnabled(bool)

    def save(self):
        """Save the input paths"""
        name = QFileDialog.getSaveFileName(None, 'Save File',
                                           ".xml")[0]
        if len(name) != 0:

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
    
            wdir.text = self.ui.wDir_lineEdit.text()
            dhm.text = self.ui.DEM_lineEdit.text()
            release.text = self.ui.release_lineEdit.text()
            infra.text = self.ui.infra_lineEdit.text()
            forest.text = self.ui.forest_lineEdit.text()
    
            tree = ET.ElementTree(root)
            tree.write(name)

    def load(self):
        xml_file = QFileDialog.getOpenFileNames(None, 'Open xml',
                                                self.directory,
                                                "xml (*.xml);;All Files (*.*)")[0]

        if len(xml_file) != 0:
            tree = ET.parse(xml_file[0])
            root = tree.getroot()
    
            try:
                self.ui.wDir_lineEdit.setText(root[0].text)
                self.directory = root[0].text
            except:
                print("No Working Directory Path in File!")
    
            try:
                self.ui.DEM_lineEdit.setText(root[1].text)
            except:
                print("No DEM Path in File!")
    
            try:
                self.ui.release_lineEdit.setText(root[2].text)
            except:
                print("No Release Path in File!")
    
            try:
                self.ui.infra_lineEdit.setText(root[3].text)
            except:
                print("No Infrastructure Path in File!")
    
            try:
                self.ui.forest_lineEdit.setText(root[4].text)
            except:
                print("No Forest Path in File!")

    def quit(self):
        QCoreApplication.quit()

    def open_wDir(self):
        """Open the Working Directory, where results are stored"""
        directory = QFileDialog.getExistingDirectory(None, 'Open Working Directory',
                                                          self.directory,
                                                          QFileDialog.ShowDirsOnly)
        if len(directory) != 0:
            self.directory = directory
            self.ui.wDir_lineEdit.setText(self.directory)

    def open_dhm(self):
        """Open digital elevation model"""
        dem_file = QFileDialog.getOpenFileNames(None, 'Open DEM',
                                                self.directory,
                                                "ascii (*.asc);;tif (*.tif);;All Files (*.*)")
        if len(dem_file[0]) != 0:
            dem = dem_file[0]
            self.ui.DEM_lineEdit.setText(dem[0])

    def open_release(self):
        """Open release layer"""
        release_file = QFileDialog.getOpenFileNames(None, 'Open Release',
                                                    self.directory,
                                                    "ascii (*.asc);;tif (*.tif);;All Files (*.*)")
        if len(release_file[0]) != 0:
            release = release_file[0]
            self.ui.release_lineEdit.setText(release[0])

    def open_infra(self):
        """Open infrastructure layer"""
        infra_file = QFileDialog.getOpenFileNames(None, 'Open Infrastructure Layer',
                                                  self.directory,
                                                  "ascii (*.asc);;tif (*.tif);;All Files (*.*)")
        if len(infra_file[0]) != 0:
            infra = infra_file[0]
            self.ui.infra_lineEdit.setText(infra[0])
            
    def open_forest(self):
        """Open forest layer"""
        forest_file = QFileDialog.getOpenFileNames(None, 'Open Forest Layer',
                                                   self.directory,
                                                   "ascii (*.asc);;tif (*.tif);;All Files (*.*)")
        forest = forest_file[0]
        self.ui.forest_lineEdit.setText(forest[0])

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
        self.calc_bool = False

        # Check if input is ok
        if self.ui.wDir_lineEdit.text() == '':
            self.showdialog('Working Directory')
            return
        if self.ui.DEM_lineEdit.text() == '':
            self.showdialog('DEM Layer')
            return
        if self.ui.release_lineEdit.text() == '':
            self.showdialog('Release Layer')
            return
        # Disable all input line Edits and Buttons
        self.set_gui_bool(False)

        # Create result directory
        time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        try:
            os.makedirs(self.ui.wDir_lineEdit.text() + '/res_{}/'.format(time_string))
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
            dem, header = io.read_raster(self.ui.DEM_lineEdit.text())
            logging.info('DEM File: {}'.format(self.ui.DEM_lineEdit.text()))
        except FileNotFoundError:
            print("Wrong filepath or filename")
            self.set_gui_bool(True)
            return

        try:
            release, release_header = io.read_raster(self.ui.release_lineEdit.text())
            logging.info('Release File: {}'.format(self.ui.release_lineEdit.text()))
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
            infra, infra_header = io.read_raster(self.ui.infra_lineEdit.text())
            if header['ncols'] == infra_header['ncols'] and header['nrows'] == infra_header['nrows']:
                print("Infra Layer ok!")
                self.calc_bool = True
                logging.info('Infrastructure File: {}'.format(self.ui.infra_lineEdit.text()))
            else:
                print("Error: Infra Layer doesn't match DEM!")
                self.set_gui_bool(True)
                return
        except:
            infra = np.zeros_like(dem)
            
        try:
            forest, forest_header = io.read_raster(self.ui.forest_lineEdit.text())
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

        alpha = self.ui.alpha_Edit.text()
        exp = self.ui.exp_Edit.text()
        flux_threshold = self.ui.flux_Edit.text()
        max_z = self.ui.z_Edit.text()
        
        logging.info('Alpha Angle: {}'.format(alpha))
        logging.info('Exponent: {}'.format(exp))
        logging.info('Flux Threshold: {}'.format(flux_threshold))
        logging.info('Max Z_delta: {}'.format(max_z))
        logging.info
        
        self.z_delta = np.zeros_like(dem)
        self.flux = np.zeros_like(dem)
        self.cell_counts = np.zeros_like(dem)
        self.z_delta_sum = np.zeros_like(dem)
        self.backcalc = np.zeros_like(dem)
        self.fp_ta = np.zeros_like(dem)
        self.sl_ta = np.zeros_like(dem)

        # Calculation
        self.calc_class = Sim.Simulation(dem, header, release, release_header, infra, forest, self.calc_bool, alpha, exp, flux_threshold, max_z)
        self.calc_class.value_changed.connect(self.update_progressBar)
        self.calc_class.finished.connect(self.thread_finished)
        logging.info('Multiprocessing starts, used cores: {}'.format(cpu_count()))
        self.calc_class.start()

    def thread_finished(self, z_delta, flux, count_array, z_delta_sum, backcalc, fp_ta, sl_ta):
        logging.info('Calculation finished, getting results.')
        for i in range(len(z_delta)):
            self.z_delta = np.maximum(self.z_delta, z_delta[i])
            self.flux = np.maximum(self.flux, flux[i])
            self.cell_counts += count_array[i]
            self.z_delta_sum += z_delta_sum[i]
            self.backcalc = np.maximum(self.backcalc, backcalc[i])
            self.fp_ta = np.maximum(self.fp_ta, fp_ta[i])
            self.sl_ta = np.maximum(self.sl_ta, sl_ta[i])
        self.output()

    def output(self):
        # time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        logging.info('Writing Output Files')
        io.output_raster(self.ui.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "flux{}".format(self.ui.outputBox.currentText()),
                         self.flux)
        io.output_raster(self.ui.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "z_delta{}".format(self.ui.outputBox.currentText()),
                         self.z_delta)
        io.output_raster(self.ui.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "FP_travel_angle{}".format(self.ui.outputBox.currentText()),
                         self.fp_ta)
        io.output_raster(self.ui.DEM_lineEdit.text(),
                         self.directory + self.res_dir + "SL_travel_angle{}".format(self.ui.outputBox.currentText()),
                         self.sl_ta)
        if not self.calc_bool:
            io.output_raster(self.ui.DEM_lineEdit.text(),
                             self.directory + self.res_dir + "cell_counts{}".format(self.ui.outputBox.currentText()),
                             self.cell_counts)
            io.output_raster(self.ui.DEM_lineEdit.text(),
                             self.directory + self.res_dir + "z_delta_sum{}".format(self.ui.outputBox.currentText()),
                             self.z_delta_sum)
        if self.calc_bool:
            io.output_raster(self.ui.DEM_lineEdit.text(),
                             self.directory + self.res_dir + "backcalculation{}".format(self.ui.outputBox.currentText()),
                             self.backcalc)

        print("Calculation finished")
        end = datetime.now().replace(microsecond=0)
        logging.info('Calculation needed: ' + str(end - self.start) + ' seconds')

        # Handle GUI
        #self.ui.progressBar.setValue(100)
        self.set_gui_bool(True)


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
        
    if 'forest' in kwargs:
        forest_path = kwargs.get('forest')
        print(forest_path)
    else:
        forest_path = None
        
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

    print("Starting...")
    print("...")

    start = datetime.now().replace(microsecond=0)
    calc_bool = False
    # Create result directory
    time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    try:
        os.makedirs(directory + '/res_{}/'.format(time_string))
        res_dir = ('/res_{}/'.format(time_string))
    except FileExistsError:
        res_dir = ('/res_{}/'.format(time_string))

    # Setup logger
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=(directory + res_dir + 'log_{}.txt').format(time_string),
                        filemode='w')

    # Start of Calculation
    logging.info('Start Calculation')
    logging.info('Alpha Angle: {}'.format(alpha))
    logging.info('Exponent: {}'.format(exp))
    logging.info('Flux Threshold: {}'.format(flux_threshold))
    logging.info('Max Z_delta: {}'.format(max_z))
    logging.info
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
            calc_bool = True
            logging.info('Infrastructure File: {}'.format(infra_path))
        else:
            print("Error: Infra Layer doesn't match DEM!")
            return
    except:
        infra = np.zeros_like(dem)
        
    try:
        forest, forest_header = io.read_raster(forest_path)
        if header['ncols'] == forest_header['ncols'] and header['nrows'] == forest_header['nrows']:
            print("Forest Layer ok!")
            prot_for_bool = True
            logging.info('Forest File: {}'.format(forest_path))
        else:
            print("Error: Forest Layer doesn't match DEM!")
            return
    except:
        forest = np.zeros_like(dem)

    logging.info('Files read in')

    z_delta = np.zeros_like(dem)
    flux = np.zeros_like(dem)
    cell_counts = np.zeros_like(dem)
    z_delta_sum = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)
    fp_ta = np.zeros_like(dem)
    sl_ta = np.zeros_like(dem)

    avaiable_memory = psutil.virtual_memory()[1]
    needed_memory = sys.getsizeof(dem)

    max_number_procces = int(avaiable_memory / (needed_memory * 10))

    print(
        "There are {} Bytes of Memory avaiable and {} Bytes needed per process. \nMax. Nr. of Processes = {}".format(
            avaiable_memory, needed_memory*10, max_number_procces))

    # Calculation
    logging.info('Multiprocessing starts, used cores: {}'.format(cpu_count()))

    if calc_bool:
        release_list = fc.split_release(release, release_header, min(mp.cpu_count() * 2, max_number_procces))
        print("Calculation 1")
        print("{} Processes started.".format(len(release_list)))
        pool = mp.Pool(len(release_list))
        results = pool.map(fc.calculation,
                           [[dem, header, infra, forest, release_pixel, alpha, exp, flux_threshold, max_z]
                            for release_pixel in release_list])
        pool.close()
        pool.join()
    else:
        release_list = fc.split_release(release, release_header, min(mp.cpu_count() * 4, max_number_procces))
        print("Calculation 2")
        print("{} Processes started.".format(len(release_list)))
        pool = mp.Pool(mp.cpu_count())
        # results = pool.map(gc.calculation, iterable)
        results = pool.map(fc.calculation_effect,
                           [[dem, header, forest, release_pixel, alpha, exp, flux_threshold, max_z] for
                            release_pixel in release_list])
        pool.close()
        pool.join()

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
    if not calc_bool:  # if no infra
        io.output_raster(dem_path,
                         directory + res_dir + "cell_counts{}".format(output_format),
                         cell_counts)
        io.output_raster(dem_path,
                         directory + res_dir + "z_delta_sum{}".format( output_format),
                         z_delta_sum)
    if calc_bool:  # if infra
        io.output_raster(dem_path,
                         directory + res_dir + "backcalculation{}".format(output_format),
                         backcalc)

    print("Calculation finished")
    print("...")
    end = datetime.now().replace(microsecond=0)
    logging.info('Calculation needed: ' + str(end - start) + ' seconds')


if __name__ == '__main__':
    #mp.set_start_method('spawn') # used in Windows
    argv = sys.argv[1:]
    #argv = ['--gui']
    #argv = ['25', '8', './examples/forest/', './examples/forest/parabola.asc', './examples/forest/release.tif', 'forest=./examples/forest/forest.tif']
    #argv = ['25', '8', './examples/forest/', './examples/forest/parabola.asc', './examples/forest/release.tif']
    if len(argv) < 1:
    	print("Too few input arguments!!!")
    	sys.exit(1)
    if len(argv) == 1 and argv[0] == '--gui':
        Flow_Py_EXEC()
    else:
        args=[arg for arg in argv if arg.find('=')<0]
        kwargs={kw[0]:kw[1] for kw in [ar.split('=') for ar in argv if ar.find('=')>0]}

        main(args, kwargs)
    
# example dam: python3 main.py 25 8 ./examples/dam/ ./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc ./examples/dam/release_dam.tif
# example dam/infra: python3 main.py 25 8 ./examples/dam/ ./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc ./examples/dam/release_dam.tif infra=./examples/dam/infra.tif flux=0.003 max_z=270
# example dam/forest: python3 main.py 25 8 ./examples/dam/ ./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc ./examples/dam/release_dam.tif forest=./examples/dam/infra.tif flux=0.003 max_z=270
