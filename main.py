import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import time

import subroutine.mainData as md
import subroutine.mainForward as mf

class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        super(MainWindow, self).__init__(parent)
        self.tab_main()
        self.information_window()

    def information_window(self):
        self.name_program = "3D Seismology Data"
        self.setWindowTitle(self.name_program)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

    def tab_main(self):
        # Initialize tab screen
        self.tabs = QTabWidget()
        self.mainData = md.WindowPick()
        self.mainForward = mf.WindowForward()

        # Add tabs
        self.tabs.addTab(self.mainData, "Generated Data")
        self.tabs.tabBar().setTabButton(0, QTabBar.RightSide, None)
        self.tabs.addTab(self.mainForward, "Forward Pseudobending")
        self.tabs.tabBar().setTabButton(0, QTabBar.RightSide, None)

        self.setCentralWidget(self.tabs)
    
    

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory().create("fusion"))
    main = MainWindow()
    main.showMaximized()
    sys.exit(app.exec_())