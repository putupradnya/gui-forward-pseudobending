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

# import mayavi
from mayavi import mlab
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor

#The actual visualization
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):

        model = np.load('D:/Project/trainee/ModelVelocity2.npy')

        mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(model),
                                    plane_orientation='x_axes',
                                    slice_index=20,
                                )
        mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(model),
                                    plane_orientation='y_axes',
                                    slice_index=20,
                                )
        mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(model),
                                    plane_orientation='z_axes',
                                    slice_index=20,
                                )
        mlab.outline()
        mlab.scalarbar(orientation='vertical')
        mlab
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some
        # VTK features require a GLContext.

        # We can do normal mlab calls on the embedded scene.
        # self.scene.mlab.test_points3d()

    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True # We need this to resize with the parent widget
                )


class WindowPick(QMainWindow):
    def __init__(self, parent=None):
        super(WindowPick, self).__init__(parent)
        self.setWindowTitle('Data Visualization')

        # Statusbar
        self.statusBar().showMessage('Nothing')

        separator = QLabel('    ')

        # Groupbox1 Create Model
        lb_x = QLabel('X Length')
        lb_y = QLabel('Y Length')
        lb_z = QLabel('Z Length')
        lb_grid = QLabel('Grid Model')
        lb_minVel = QLabel('Min Velocity')
        lb_maxVel = QLabel('Max Velocity')
        self.bt_outModel = QPushButton('Load')
        self.bt_outModel.clicked.connect(self.saveModel)
        self.bt_processModel = QPushButton('Create Model')
        self.bt_processModel.clicked.connect(self.createModel)
        self.bt_processModel.setFixedSize(100, 25)

        self.lb_x = QLineEdit()
        self.lb_y = QLineEdit()
        self.lb_z = QLineEdit()
        self.lb_grid = QLineEdit()
        self.minVel = QLineEdit()
        self.maxVel = QLineEdit()
        self.lb_outModel = QLineEdit()
        self.lb_outModel.setPlaceholderText('Click load to save model (*.npy)')

        modelform = QFormLayout()
        modelform.addRow(lb_x, self.lb_x)
        modelform.addRow(lb_y, self.lb_y)
        modelform.addRow(lb_z, self.lb_z)
        modelform.addRow(lb_grid, self.lb_grid)
        modelform.addRow(lb_minVel, self.minVel)
        modelform.addRow(lb_maxVel, self.maxVel)
        modelform.addRow(self.bt_outModel, self.lb_outModel)
        modelform.addRow(separator, self.bt_processModel)

        
        grupCreateModel = QGroupBox('Create Model')
        grupCreateModel.setFixedWidth(300)
        grupCreateModel.setLayout(modelform)

        # Grupbox2 Create Gaussian Anomaly
        bt_loadModel = QPushButton('Load')
        lb_min_x = QLabel('Min X')
        lb_min_y = QLabel('Min Y')
        lb_min_z = QLabel('Min Z')
        lb_max_x = QLabel('Min X')
        lb_max_y = QLabel('Min Y')
        lb_max_z = QLabel('Min Z')
        lb_value_anomaly = QLabel('Anomaly Value')

        self.loadModel = QLineEdit()
        self.lb_min_x = QLineEdit()
        self.lb_min_y = QLineEdit()
        self.lb_min_z = QLineEdit()
        self.lb_max_x = QLineEdit()
        self.lb_max_y = QLineEdit()
        self.lb_max_z = QLineEdit()
        self.lb_value_anomaly = QLineEdit()

        bt_add = QPushButton('Add Anomaly')
        bt_add.setFixedSize(100, 25)

        anomalyform = QFormLayout()
        anomalyform.addRow(bt_loadModel, self.loadModel)
        anomalyform.addRow(lb_min_x, self.lb_min_x)
        anomalyform.addRow(lb_min_y, self.lb_min_y)
        anomalyform.addRow(lb_min_z, self.lb_min_z)
        anomalyform.addRow(lb_max_x, self.lb_max_x)
        anomalyform.addRow(lb_max_y, self.lb_max_y)
        anomalyform.addRow(lb_max_z, self.lb_max_z)
        anomalyform.addRow(lb_value_anomaly, self.lb_value_anomaly)
        anomalyform.addRow(separator, bt_add)

        grupAnomali = QGroupBox('Create Gaussian Anomaly')
        grupAnomali.setFixedWidth(300)
        grupAnomali.setLayout(anomalyform)


        # Grupbox3 Display model
        grupDisplay= QGroupBox('Display Model')
        grupDisplay.setFixedWidth(300)
        
        # Layout toolbox
        vbox = QVBoxLayout()
        vbox.addWidget(grupCreateModel)
        vbox.addWidget(grupAnomali)
        vbox.addWidget(grupDisplay)

        # Connect to Mayavi
        self.visualmayavi = Visualization()
        self.disp_mayavi = self.visualmayavi.edit_traits().control
        self.disp_mayavi.setParent(self)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.disp_mayavi)

        # control layout
        ui = QHBoxLayout()
        ui.addLayout(vbox)
        ui.addLayout(layout)

        self.widget = QWidget()
        self.widget.setLayout(ui)
        self.setCentralWidget(self.widget)

    def saveModel(self):
        self.pathOutData = QFileDialog.getExistingDirectory(self, 'Select directory to save model')
        self.lb_outModel.setText(str(self.pathOutData))

    def createModel(self):
        xmax = int(self.lb_x.text())
        ymax = int(self.lb_y.text())
        zmax = int(self.lb_z.text())
        grid = int(self.lb_grid.text())
        vmin = np.float(self.minVel.text())
        vmax = np.float(self.maxVel.text())

        # Initial model
        x = np.r_[0:xmax:grid]
        y = np.r_[0:ymax:grid]
        z = np.r_[0:zmax:grid]

        v = np.linspace(vmin, vmax, len(z))

        Model = np.zeros((len(x), len(y), len(z)))
        for i in range(len(x)):
            for j in range(len(y)):
                for k in range(len(z)):
                    Model[i][j][k] = v[k]

        # Create output data
        path = os.path.join(str(self.pathOutData) + '/Velocity Model')
        np.save(os.path.join(path), Model)

        self.statusBar().showMessage('Successfully create model')
        QMessageBox.about(self, 'Information', 'Model created successfully!')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = WindowPick()
    w.showMaximized()
    sys.exit(app.exec_())