# from relativePath import resourcePath
import os
import sys

def resourcePath(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        path = os.getcwd()
        base_path = os.path.join(path,'addins')
    return os.path.join(base_path,relative_path)

# tab icon
icon_window = resourcePath('SeiraAppBut.png')
icon_pick = resourcePath('picking.png')
icon_depth = resourcePath('depth.png')
icon_tomo = resourcePath('tomo.png')
icon_jointomo = resourcePath('joint.png')
icon_simulation = resourcePath('icons8-solidworks-flow-simulation-100.png')
icon_load = resourcePath('open.png')
icon_save = resourcePath('save.png')
icon_import = resourcePath('import.png')
icon_pick2 = resourcePath('pick.png')
imgSplashScreen = resourcePath('seiraApp.png')