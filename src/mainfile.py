import sys
import os
import csv
from PyQt5 import QtGui
from numpy import polyval, poly1d, polyfit
import numpy as np
import pandas as pd
import pyqtgraph as pg
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
import scipy.signal
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QFileDialog, QLineEdit, QGraphicsScene, QLabel, \
    QComboBox, QColorDialog
from PyQt5.QtGui import QIcon, QPixmap
from pathlib import Path
import math
import pyqtgraph.exporters
from PyQt5.QtWidgets import QBoxLayout, QLineEdit, QSpinBox
from PyQt5.QtWidgets import QLabel
from PyQt5.QtGui import QDoubleValidator, QValidator
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from scipy.interpolate import interp1d
import numpy.fft as fft
from scipy.interpolate import make_interp_spline
from numpy import savetxt
import scipy
from scipy import interpolate
from design import Ui_MainWindow, MplCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib as mpl
import statistics
from sympy import S, symbols, printing

# from IPython.display import display, Latex


class MainWindow(qtw.QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.AmplitudeList = []
        self.TimeList = []
        # self.array_of_chunksx=[]
        # self.array_of_chunksy=[]
        self.number_of_chunk = 0
        #self.playTimer = QtCore.QTimer()
        self.pen1 = pg.mkPen(color=(255, 0, 0))
        self.pen2 = pg.mkPen(color=(0, 139, 139))
        self.pen3 = pg.mkPen(color=(255, 255, 0))
        self.pen4 = pg.mkPen(color=(0, 255, 0))
        self.pen5 = pg.mkPen(color=(138, 43, 226))
        self.pen6 = pg.mkPen(color=(255, 99, 71))
        self.pen7 = pg.mkPen(color=(128, 128, 128))
        self.pen8 = pg.mkPen(color=(0, 128, 128))
        self.pen9 = pg.mkPen(color=(255, 20, 147))
        self.pen10 = pg.mkPen(color=(118, 238, 198))
        self.Color = [self.pen1, self.pen2, self.pen3, self.pen4,
                      self.pen5, self.pen6, self.pen7, self.pen8, self.pen9, self.pen10]
        self.order_interpolation_val = 0
        self.ui.actionOpen.triggered.connect(lambda: self.open_file())

        self.ui.mapping_pushButton.clicked.connect(lambda: self.show_mapping())
        self.ui.startmap_pushButton.clicked.connect(
            lambda: self.error_mapping())
        self.ui.cancelmap_pushButton.clicked.connect(
            lambda: self.cancel_mapping())
        self.ui.pushButton_extrapolate.clicked.connect(
            lambda: self.extrapolate())
        self.ui.addEquation_pushButton.clicked.connect(lambda: self.addequ())
        self.ui.showeqn_pushButton.clicked.connect(lambda: self.write())
        self.ui.spinBox_fitchunk.valueChanged.connect(
            self.data_in_multichuncks)

        self.ui.spinBox_mapchunk.hide()
        self.ui.spinBox_mappoly.hide()
        self.ui.spinBox_mapoverlap.hide()
        self.ui.startmap_pushButton.hide()
        self.ui.cancelmap_pushButton.hide()
        self.ui.xaxis_comboBox.hide()
        self.ui.yaxis_comboBox.hide()
        self.ui.chunk_no_map_label.hide()
        self.ui.poly_no_map_label.hide()
        self.ui.overlap_map_label.hide()
        self.ui.xaxis_label.hide()
        self.ui.yaxis_label.hide()
        self.ui.progressBar.hide()
        self.ui.line_5.hide()
        self.equations1 = []
        self.equations2 = []
        self.equations3 = []
        self.equations4 = []
        self.coffs = []
        self.c1 = ""
        # self.c2=[]

    def open_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.file_url, _ = QFileDialog.getOpenFileName(None, "QFileDialog.getOpenFileName()", "",
                                                       "All Files (*);;csv Files (*.csv)", options=options)
        if self.file_url:
            print(self.file_url)
            filename = Path(self.file_url).stem
            print(filename)
            self.read_file(self.file_url)

    def read_file(self, file_path):
        self.ui.fitting_graph.clear()
        path = file_path

        read_data = pd.read_csv(path)

        self.amplitude = read_data.values[:, 1]
        self.AmplitudeList.append(self.amplitude)
        self.time = read_data.values[:, 0]
        self.TimeList.append(self.time)
        for i in range(len(self.AmplitudeList)):
            if i <= 6:
                self.ui.fitting_graph.plot(
                    self.TimeList[i], self.AmplitudeList[i], pen=self.Color[5])

    def get_number_of_chunk(self):
        self.number_of_chunk = self.ui.spinBox_fitchunk.value()

        return (self.number_of_chunk)

    def get_order_interpolation(self):
        self.order_interpolation_val = self.ui.spinBox_fitorder.value()
        return (self.order_interpolation_val)

    def data_in_multichuncks(self):

        self.array_of_chunksx = np.array(
            np.array_split(self.time, self.get_number_of_chunk()))
        self.array_of_chunksy = np.array(np.array_split(
            self.amplitude, self.get_number_of_chunk()))
        for i in range(self.get_number_of_chunk()):
            self.ui.fitting_graph.plot(
                self.array_of_chunksx[i], self.array_of_chunksy[i], pen=self.Color[i])

    def write_equation(self, mathTex):
        fig = mpl.figure.Figure()
        fig.patch.set_facecolor('none')
        fig.set_canvas(FigureCanvasAgg(fig))
        renderer = fig.canvas.get_renderer()
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        ax.patch.set_facecolor('none')
        t = ax.text(0, 0, r"$%s$" % (mathTex),
                    ha='left', va='bottom', fontsize=10)

    # ---- fit figure size to text artist ----

        fwidth, fheight = fig.get_size_inches()
        fig_bbox = fig.get_window_extent(renderer)
        text_bbox = t.get_window_extent(renderer)
        tight_fwidth = text_bbox.width * fwidth / fig_bbox.width
        tight_fheight = text_bbox.height * fheight / fig_bbox.height
        fig.set_size_inches(tight_fwidth, tight_fheight)

# ---- convert mpl figure to QPixmap ----

        buf, size = fig.canvas.print_to_buffer()
        qimage = QtGui.QImage.rgbSwapped(QtGui.QImage(buf, size[0], size[1],
                                                      QtGui.QImage.Format_ARGB32))
        qpixmap = QtGui.QPixmap(qimage)
        return qpixmap

    def order_interpolation(self):
        self.ui.fitting_graph.clear()
        self.read_file(self.file_url)
        self.data_in_multichuncks()
        for i in range(self.get_number_of_chunk()):
            self.coeff1 = np.polyfit(
                self.array_of_chunksx[i],  self.array_of_chunksy[i], self.get_order_interpolation())
            print(self.coeff1)
            self.coffs.append(self.coeff1)
            self.model1 = np.poly1d(np.polyfit(
                self.array_of_chunksx[i],  self.array_of_chunksy[i], self.get_order_interpolation()))
            self.point = self.array_of_chunksx[i]
            self.end = self.point[-1]
            self.start = self.point[0]
            self.polyline = np.linspace(
                self.start, self.end, len(self.array_of_chunksx[i]))
            self.ui.fitting_graph.plot(self.polyline, self.model1(
                self.polyline), symbol='o', pen=self.Color[i])
            self.error1 = self.sqr_error(
                self.model1, self.array_of_chunksx[i], self.array_of_chunksy[i],)
        print("ccccccccccccccccccc", self.coffs)

    def write(self):
        self.order_interpolation()
        self.coefficients = self.coffs[self.ui.comboBox_chunk_eqn.currentIndex(
        )-1]
        self.c1 = "chunk{}:".format(
            self.ui.comboBox_chunk_eqn.currentIndex()+1)
        for i in range(len(self.coefficients)-1):
            coeff = str(round(self.coefficients[i], 2))
            if coeff[0] == "-":
                self.c1 += r"{}$x$^{}".format(
                    round(self.coefficients[i], 2), len(self.coefficients)-1-i)
            else:
                self.c1 += r"+{}$x$^{}".format(
                    round(self.coefficients[i], 2), len(self.coefficients)-1-i)
        lastcoeff = str(round(self.coefficients[-1], 2))
        if lastcoeff[0] == "-":
            self.c1 += r"{}".format(lastcoeff)
        else:
            self.c1 += r"+ {}".format(lastcoeff)
        self.ui.showfittingequation_label.setPixmap(
            self.write_equation(self.c1))
        self.ui.showerror_label.setPixmap(
            self.write_equation(r"{}".format(round(self.error1, 3))))
        print(self.c1)

    def addequ(self):
        for i in range(self.get_number_of_chunk()):
            self.ui.comboBox_chunk_eqn.addItem("chunk {}".format(i+1))

    def sqr_error(self, p, xi, yi):

        self.diff2 = (p(xi)-yi)**2

        return self.diff2.sum()

    def show_mapping(self):
        self.ui.spinBox_mapchunk.show()
        self.ui.spinBox_mapoverlap.show()
        self.ui.spinBox_mappoly.show()
        self.ui.startmap_pushButton.show()
        self.ui.xaxis_comboBox.show()
        self.ui.yaxis_comboBox.show()
        self.ui.chunk_no_map_label.show()
        self.ui.poly_no_map_label.show()
        self.ui.overlap_map_label.show()
        self.ui.xaxis_label.show()
        self.ui.yaxis_label.show()
        self.ui.line_5.show()
        self.ui.progressBar.show()
        self.ui.progressBar.setValue(0)

    def cancel_mapping(self):
        self.ui.progressBar.setValue(0)
        # LOOP TO DELETE THE OLD WIDGET>>SPECTROGRAM<< WITH ITS ALL ITEMS
        for i in reversed(range(self.ui.mapping_layout.count())):
            self.ui.mapping_layout.itemAt(i).widget().deleteLater()
        self.ui.cancelmap_pushButton.hide()
        self.ui.startmap_pushButton.show()

    def error_mapping(self):
        # make an array from 0 to ovalap value
        self.overlap_array = np.arange(
            0, self.ui.spinBox_mapoverlap.value(), 1)
        # make an array from 0 to poly value
        self.polynomial_array = np.arange(
            0, self.ui.spinBox_mappoly.value(), 1)
        # make an array from 0 to chunk n value
        self.chunk_n_array = np.arange(0, self.ui.spinBox_mapchunk.value(), 1)

        self.ui.startmap_pushButton.hide()
        self.ui.cancelmap_pushButton.show()
        # LOOP TO DELETE THE OLD WIDGET>>SPECTROGRAM<< WITH ITS ALL ITEMS
        for i in reversed(range(self.ui.mapping_layout.count())):
            self.ui.mapping_layout.itemAt(i).widget().deleteLater()

        if (self.ui.xaxis_comboBox.currentIndex() == 0 and self.ui.yaxis_comboBox.currentIndex() == 0):  # poly & overlap

            self.error_array = [np.zeros(
                self.ui.spinBox_mappoly.value())]*int(self.ui.spinBox_mapoverlap.value())
            self.completed = 0      # % of progress bar
            step = int(100/(self.ui.spinBox_mappoly.value() *
                       self.ui.spinBox_mapoverlap.value()))  # step of progress bar
            self.ui.progressBar.show()
            for j in range(self.ui.spinBox_mapoverlap.value()):  # loop for rows
                for i in range(self.ui.spinBox_mappoly.value()):  # loop for colums
                    # make the fitting process for each colum value(poly number)
                    p = polyfit(self.time, self.amplitude, int(i))
                    # the y array after fitted
                    y_fitting = polyval(p, self.time)
                    # modify-------------------------------------------------
                    error = abs((statistics.stdev(
                        self.amplitude)-statistics.stdev(y_fitting))/statistics.stdev(self.amplitude))*100
                    self.error_array[j][i] = error
                    print(str(error))
                    # self.ui.tableWidget.setItem(j,i,QtWidgets.QTableWidgetItem(str(error)))#print avg(the error val) in table
                    self.completed += step  # increase the progress bar
                    self.ui.progressBar.setValue(self.completed)
                if(j == self.ui.spinBox_mapoverlap.value()-1):
                    # set progress bar at 100% at final value
                    self.ui.progressBar.setValue(100)

                    # LOOP TO DELETE THE OLD WIDGET>>SPECTROGRAM<< WITH ITS ALL ITEMS
                    for i in reversed(range(self.ui.mapping_layout.count())):
                        print(i)
                        if i == 2:
                            self.ui.mapping_layout.itemAt(
                                i).widget().deleteLater()
                        elif i == 1:
                            self.ui.mapping_layout.itemAt(
                                i).widget().deleteLater()

                    self.map_canvas = MplCanvas(
                        self.ui.centralwidget, width=10, height=10, dpi=100)
                    self.ui.mapping_layout.addWidget(self.map_canvas)
                    self.map_canvas.axes.imshow(
                        self.error_array, cmap='jet', aspect='auto', origin="lower")
                    im = self.map_canvas.axes.imshow(
                        self.error_array, cmap='jet', aspect='auto', origin="lower")
                    self.map_canvas.axes.set_xlabel("Polynomial no.")
                    self.map_canvas.axes.xaxis.set_label_coords(1.17, -0.035)
                    self.map_canvas.axes.set_ylabel("Overlapping")
                    self.map_canvas.figure.colorbar(
                        im, label="Error %", orientation="vertical")
                    self.map_canvas.draw()

        elif (self.ui.xaxis_comboBox.currentIndex() == 0 and self.ui.yaxis_comboBox.currentIndex() == 1):  # poly & chunks number

            self.error_array = [np.zeros(
                self.ui.spinBox_mappoly.value())]*int(self.ui.spinBox_mapchunk.value())
            self.completed = 0
            step = int(100/(self.ui.spinBox_mappoly.value()
                       * self.ui.spinBox_mapchunk.value()))
            self.ui.progressBar.show()
            for j in range(1, self.ui.spinBox_mapchunk.value()+1):
                self.start = 0
                self.end = 0
                error_array_list = []
                for c in range(1, j):
                    self.chunck_size = int(len(self.time)/j)
                    self.start = self.end
                    self.end += self.chunck_size
                    self.rangex = self.time[self.start:self.end:1]
                    self.rangey = self.amplitude[self.start:self.end:1]
                    error_array = []
                    for i in range(1, self.ui.spinBox_mappoly.value()+1):
                        p = polyfit(self.rangex, self.rangey, int(i))
                        f = polyval(p, self.rangex)
                        error = abs((statistics.stdev(
                            self.rangey)-statistics.stdev(f))/statistics.stdev(self.rangey))*100
                        error_array.append(error)
                    error_array_list.append(error_array)
                list_of_chunk_parts = []
                for part in error_array_list:
                    list_of_chunk_parts = list_of_chunk_parts+part
                # print(""+str(j)+str(list_of_chunk_parts))
                for Pn in range(0, j-1):
                    eachpoly_part = []
                    for Pnn in range(0, j-1):
                        eachpoly_part.append(list_of_chunk_parts[Pnn])
                    print(str(eachpoly_part))
                    self.error_array[j-1][i-1] = max(eachpoly_part)
                    self.completed += step
                    self.ui.progressBar.setValue(self.completed)
                if(j == self.ui.spinBox_mapchunk.value()):
                    print(str(self.error_array))
                    self.ui.progressBar.setValue(100)
                    # LOOP TO DELETE THE OLD WIDGET>>SPECTROGRAM<< WITH ITS ALL ITEMS
                    for i in reversed(range(self.ui.mapping_layout.count())):
                        print(i)
                        if i == 2:
                            self.ui.mapping_layout.itemAt(
                                i).widget().deleteLater()

                        elif i == 1:
                            self.ui.mapping_layout.itemAt(
                                i).widget().deleteLater()

                    self.map_canvas = MplCanvas(
                        self.ui.centralwidget, width=10, height=10, dpi=100)
                    self.ui.mapping_layout.addWidget(self.map_canvas)
                    self.map_canvas.axes.imshow(
                        self.error_array, cmap='jet', aspect='auto', origin="lower")
                    im = self.map_canvas.axes.imshow(
                        self.error_array, cmap='jet', aspect='auto', origin="lower")
                    self.map_canvas.axes.set_xlabel("Polynomial no.")
                    self.map_canvas.axes.xaxis.set_label_coords(1.13, -0.035)
                    self.map_canvas.axes.set_ylabel("Chunks no.")
                    self.map_canvas.figure.colorbar(
                        im, label="Error %", orientation="vertical")
                    self.map_canvas.draw()

        elif (self.ui.xaxis_comboBox.currentIndex() == 1 and self.ui.yaxis_comboBox.currentIndex() == 1):  # overlap & chunks number
            self.error_array = [np.zeros(
                self.ui.spinBox_mapoverlap.value()+1)]*int(self.ui.spinBox_mapchunk.value())
            self.completed = 0
            step = int(100/(self.ui.spinBox_mapoverlap.value()
                       * self.ui.spinBox_mapchunk.value()))
            self.ui.progressBar.show()
            for j in range(1, self.ui.spinBox_mapchunk.value()+1):
                if j == 1:
                    p = polyfit(self.time, self.amplitude, 5)
                    f = polyval(p, self.time)
                    error = abs((statistics.stdev(
                        self.amplitude)-statistics.stdev(f))/statistics.stdev(self.amplitude))*100
                    for i in range(self.ui.spinBox_mapoverlap.value()+1):
                        self.error_array[0][i] = error
                else:
                    start = 0
                    end = 0
                    chunck_size = int(len(self.time)/j)
                    start = end
                    end += chunck_size
                    # x points of the chunk
                    self.rangex = self.time[start:end:1]
                    # y points ............
                    self.rangey = self.amplitude[start:end:1]
                    p = polyfit(self.rangex, self.rangey, 5)
                    f = polyval(p, self.rangex)
                    error_part1 = abs((statistics.stdev(
                        self.rangey)-statistics.stdev(f))/statistics.stdev(self.rangey))*100
                    part_at_different_overlaps = [
                        np.zeros(self.ui.spinBox_mapoverlap.value())]*int(j-1)
                    endpart = end
                    for part in range(2, j):
                        start = end
                        end += chunck_size
                        for i in range(self.ui.spinBox_mapoverlap.value()+1):
                            if start == end:
                                break
                            start = start-int((i*chunck_size)/100)
                            # x points of the chunk
                            self.rangex = self.time[start:end:1]
                            # y points ............
                            self.rangey = self.amplitude[start:end:1]
                            p = polyfit(self.rangex, self.rangey, 5)
                            f = polyval(p, self.rangex)
                            error = abs((statistics.stdev(
                                self.rangey)-statistics.stdev(f))/statistics.stdev(self.rangey))*100
                            part_at_different_overlaps[int(
                                part-2)][i-1] = error
                        endpart += chunck_size
                        for i in range(self.ui.spinBox_mapoverlap.value()+1):
                            chunks_error_array = []
                            chunks_error_array.append(error_part1)
                            for part_1 in range(2, j):
                                chunks_error_array.append(
                                    part_at_different_overlaps[part_1-2][i-1])
                            self.error_array[j-1][i] = max(chunks_error_array)

                if(j == self.ui.spinBox_mapchunk.value()):
                    print(str(self.error_array))
                    self.ui.progressBar.setValue(100)
                    # LOOP TO DELETE THE OLD WIDGET>>SPECTROGRAM<< WITH ITS ALL ITEMS
                    for i in reversed(range(self.ui.mapping_layout.count())):
                        print(i)
                        if i == 2:
                            self.ui.mapping_layout.itemAt(
                                i).widget().deleteLater()
                        elif i == 1:
                            self.ui.mapping_layout.itemAt(
                                i).widget().deleteLater()

                    self.map_canvas = MplCanvas(
                        self.ui.centralwidget, width=10, height=10, dpi=100)
                    self.ui.mapping_layout.addWidget(self.map_canvas)
                    self.map_canvas.axes.imshow(
                        self.error_array, cmap='summer', aspect='auto', origin="lower")
                    im = self.map_canvas.axes.imshow(
                        self.error_array, cmap='summer', aspect='auto', origin="lower")
                    self.map_canvas.axes.set_xlabel("Overlapping")
                    self.map_canvas.axes.xaxis.set_label_coords(1.13, -0.035)
                    self.map_canvas.axes.set_ylabel("Chunks no.")
                    self.map_canvas.figure.colorbar(
                        im, label="Error %", orientation="vertical")
                    self.map_canvas.draw()

    def extrapolate(self):
        self.ui.fitting_graph.clear()
        self.ui.fitting_graph.plot(
            self.time, self.amplitude, pen=self.Color[5])

        percenatge = self.ui.spinBox_signalpercentage.value()
        order = self.ui.spinBox_fitorder.value()
        self.newamp = self.amplitude[:int(
            (percenatge/100)*len(self.amplitude))]
        self.newtime = self.time[:int((percenatge/100)*len(self.time))]
        newfit = np.polyfit(self.newtime, self.newamp, order)
        self.newval = np.polyval(newfit, self.newtime)
        self.ui.fitting_graph.plot(
            self.newtime, self.newval, symbol='o', pen=self.pen3)

        # newfit2 = np.polyfit(self.time, self.amplitude, order)
        self.newval2 = np.polyval(newfit, self.time)
        self.ui.fitting_graph.plot(
            self.time, self.newval2, symbol='o', pen=self.pen4)


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())
