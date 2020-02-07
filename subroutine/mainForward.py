import time
import os
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pathlib import Path
import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# import subroutine
import subroutine.module.readdata as rd
import subroutine.module.forward as fw
from subroutine.thread.threading import Worker, MessageBox
from subroutine.time.tictac import tic, tac

class WindowForward(QMainWindow):
    def __init__(self, parent=None):
        super(WindowForward, self).__init__(parent)
        self.setWindowTitle('Data Visualization')

        # Threading core
        self.threadpool = QThreadPool()

        # Status bar
        self.statusBar().showMessage('Nothing')

        # a figure instance to plot on
        self.figure   = plt.figure()
        self.canvas   = FigureCanvas(self.figure)
        self.toolbar  = NavigationToolbar(self.canvas, self)

        separator = QLabel('    ')

        # Groupbox input
        grupInput = QGroupBox('Forward Modeling')
        grupInput.setFixedWidth(300)

        self.bt_model = QPushButton('Load')
        self.bt_model.clicked.connect(self.openModel)

        self.bt_event = QPushButton('Load')
        self.bt_event.clicked.connect(self.openEvent)

        self.bt_station = QPushButton('Load')
        self.bt_station.clicked.connect(self.openStation)

        self.bt_dataout = QPushButton('Load')
        self.bt_dataout.clicked.connect(self.outputPath)
        
        self.bt_createmodel = QPushButton('Create Model')
        self.bt_createmodel.clicked.connect(self.createModel)
        self.bt_execute = QPushButton('EXEC!')
        self.bt_execute.clicked.connect(self.clickButtonExecute)

        lb_data     = QLabel('Input data')
        lb_param    = QLabel('Parameter')
        lb_outdata  = QLabel('Output data')
        lb_model    = QLabel('Build model')

        lb_nx               = QLabel('Length X')
        lb_ny               = QLabel('Length Y')
        lb_nz               = QLabel('Length Z')
        lb_grid             = QLabel('Grid size')
        lb_niter            = QLabel('Iteration')
        lb_cacah            = QLabel('Segment')
        lb_resolution       = QLabel('Split Resolution')
        self.lb_path_model  = QLineEdit()
        self.lb_path_model.setPlaceholderText('path of velocity model')
        self.lb_path_event  = QLineEdit()
        self.lb_path_event.setPlaceholderText('path of event')
        self.lb_path_station = QLineEdit()
        self.lb_path_station.setPlaceholderText('path of station')
        self.lb_path_outdata = QLineEdit()
        self.lb_path_outdata.setPlaceholderText('path of output data')

        self.niter          = QLineEdit()
        self.cacah          = QLineEdit()
        self.nx             = QLineEdit()
        self.ny             = QLineEdit()
        self.nz             = QLineEdit()
        self.gridsize       = QLineEdit()
        self.split_resolution = QLineEdit()

        fboxMeth = QFormLayout()
        fboxMeth.addRow(lb_model)
        fboxMeth.addRow(lb_nx, self.nx)
        fboxMeth.addRow(lb_ny, self.ny)
        fboxMeth.addRow(lb_nz, self.nz)
        fboxMeth.addRow(lb_grid, self.gridsize)
        fboxMeth.addRow(separator)

        fboxMeth.addRow(lb_data)
        fboxMeth.addRow(self.bt_model, self.lb_path_model)
        fboxMeth.addRow(self.bt_event, self.lb_path_event)
        fboxMeth.addRow(self.bt_station, self.lb_path_station)
        fboxMeth.addRow(separator, self.bt_createmodel)

        fboxMeth.addRow(lb_param)
        fboxMeth.addRow(lb_niter, self.niter)
        fboxMeth.addRow(lb_cacah, self.cacah)
        fboxMeth.addRow(lb_resolution, self.split_resolution)
        fboxMeth.addRow(self.bt_dataout, self.lb_path_outdata)
        fboxMeth.addRow(separator, self.bt_execute)
        grupInput.setLayout(fboxMeth)

        grupDisplay = QGroupBox('Display')
        
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        grupDisplay.setLayout(layout)

        # Tab layout
        hbox = QHBoxLayout(self)
        hbox.addWidget(grupInput)
        hbox.addWidget(grupDisplay)

        self.widget = QWidget()
        self.widget.setLayout(hbox)
        # self.widget.setLayout(self.statbar)
        self.setCentralWidget(self.widget)
    
    def outputPath(self):
      self.pathOutData = QFileDialog.getExistingDirectory(self, 'Select directory to save data')
      self.lb_path_outdata.setText(str(self.pathOutData))
    
    def openEvent(self):
      path = QFileDialog.getOpenFileName(self, 'Open event data', filter='*.event')
      path_event = Path(path[0])
      self.lb_path_event.setText(str(path_event))
      
      self.xevent, self.yevent, self.zevent = rd.readEvent(path_event)
      
      # Status bar
      self.statusBar().showMessage(str(path_event))
    
    def openStation(self):
      path = QFileDialog.getOpenFileName(self, 'Open event data', filter='*.stat')
      path_station = Path(path[0])
      self.lb_path_station.setText(str(path_station))

      self.xstation, self.ystation, self.zstation = rd.readStation(path_station)

    
    def openModel(self):
      path = QFileDialog.getOpenFileName(self, 'Open event data', filter='*.npy')
      path_model = Path(path[0])
      self.lb_path_model.setText(str(path_model))

      self.velocityModel = np.load(str(path_model))

    def createModel(self):
      # create an axis
      self.figure.clear()
      self.ax = self.figure.add_subplot(111, projection='3d')
      self.x = np.r_[0 : int(self.nx.text()) : np.float(self.gridsize.text())]
      self.y = np.r_[0 : int(self.ny.text()) : np.float(self.gridsize.text())]
      self.z = np.r_[0 : int(self.nz.text()) : np.float(self.gridsize.text())]

      self.ax.scatter(self.xevent, self.yevent, self.zevent, color='red',marker='o', linewidth=0.1, label='Event')
      self.ax.scatter(self.xstation, self.ystation, self.zstation, color='blue', marker='v', linewidth=0.1, label='Station')
      plt.legend()
      self.ax.invert_zaxis()
      plt.ion()

      self.figure.subplots_adjust(left=0.037,right=0.988,
                      bottom=0.067,top=0.9,
                      hspace=0.2,wspace=0.2)
      self.ax.invert_yaxis()
      self.ax.set_title('Data Visualization')
      self.ax.set_ylabel('Y Coordinate')
      self.ax.set_xlabel('X Coordinate')
      self.ax.set_zlabel('Depth (m)')
      self.ax.set_xlim(np.max(self.x))
      self.ax.set_xlim(np.max(self.y))
      self.ax.set_xlim(np.max(self.z))
      plt.tight_layout(pad = 1.5, w_pad=1.5, h_pad=1)
      QMessageBox.about(self, 'Information', 'Data loaded successfully!')

    def running(self, progress_callback):
      progress_callback.emit('Running data')
      niter = int(self.niter.text())
      cacah = int(self.cacah.text())
      split_resolution = int(self.split_resolution.text())

      for i in range(len(self.xevent)):
        xs = self.xevent[i]
        ys = self.yevent[i]
        zs = self.zevent[i]

        for j in range(len(self.xstation)):
          xr = self.xstation[j]
          yr = self.ystation[j]
          
          xb, yb, zb, time_total = fw.forwardProcess(self, niter, cacah, split_resolution, self.x, self.y, self.z, self.velocityModel,xs,ys,zs,xr,yr)
          self.statusBar().showMessage('Running Event: ' + str(i) + ' of ' + str(j))
          print('Running Event: ' + str(i) + ' of ' + str(len(self.xevent)) 
                + ' in station : ' + str(j) + ' from : ' + str(len(self.xstation)) + ' stations'
                + ' time elapsed: ' + str(time_total) + ' sec.')

    def executeClickButtonExecute(self, progress_callback):
      self.statusBar().showMessage('Processing data')
      tic()
      
      # self.running(progress_callback)
      niter = int(self.niter.text())
      cacah = int(self.cacah.text())
      split_resolution = int(self.split_resolution.text())
      time_ray = []

      for i in range(len(self.xevent)):
        xs = self.xevent[i]
        ys = self.yevent[i]
        zs = self.zevent[i]

        for j in range(len(self.xstation)):
          xr = self.xstation[j]
          yr = self.ystation[j]
          
          xb, yb, zb, time_total = fw.forwardProcess(self, niter, cacah, split_resolution, self.x, self.y, self.z, self.velocityModel,xs,ys,zs,xr,yr)
          self.ax.plot(xb, yb, zb, linewidth=0.1, color='b', label='Raypath')

          plt.ion()
          time_ray.append(time_total)
          self.statusBar().showMessage('Running Event: ' + str(i) + ' of ' + str(len(self.xevent)) 
                + ' in station : ' + str(j) + ' from  ' + str(len(self.xstation)) + ' stations'
                + ' || time elapsed: ' + str(time_total) + ' Sec')
      
      elt = tac()

      # Create output data
      path = os.path.join(str(self.pathOutData) + '/result.forward')
      file = open(os.path.join(path), 'w')
      file.write('No \t' + 'XS \t' + 'YS \t' + 'ZS\t' 
                  + 'XS \t' + 'YS\t' + 'ZS \t' + 'Time(s) \n')
      no = 0
      print(' error line 242')
      for i in range(len(self.xevent)):
        xs = self.xevent[i]
        ys = self.yevent[i]
        zs = self.zevent[i]

        for j in range(len(self.xstation)):
          xr = self.xstation[j]
          yr = self.ystation[j]
          zr = 0

          file.write(str(no) +'\t' + str(xs) + '\t' + str(ys) + '\t' + str(zs)
                      + str(xr) + '\t' + str(yr) + '\t' + str(zr) + '\t' + '%0.4f' %(time_ray[no]) + '\n')
        
          no+=1
      
      file.close()

      return elt

    def outputClickButtonExecute(self, s):
        # file = open('time elapsed', 'w')
        # file.write('Elapsed Time:' + '\t' + s + '\n')
      self.statusBar().showMessage('Finish')

    def errorClickButtonExecute(self):
      self.message = MessageBox()
      self.message.show()

    def completeClickButtonExecute(self):
      self.statusBar().showMessage('Finish')
        # print('Status: Finished')

    def progressClickButtonExecute(self, n):
      self.statusBar().showMessage(n)
        # print('error')
    def clickButtonExecute(self):
      # Pass the function to execute
      worker = Worker(
          self.executeClickButtonExecute)  # Any other args, kwargs are passed to the run function
      worker.signals.result.connect(self.outputClickButtonExecute)
      worker.signals.error.connect(self.errorClickButtonExecute)
      worker.signals.finished.connect(self.completeClickButtonExecute)
      worker.signals.progress.connect(self.progressClickButtonExecute)

      # Execute
      self.threadpool.start(worker)
      
if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = WindowForward()
    w.showMaximized()
    sys.exit(app.exec_())