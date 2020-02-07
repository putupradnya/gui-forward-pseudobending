import sys
import os

def resourcePath(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        path = os.getcwd()
        base_path = os.path.join(path,'addins')
    return os.path.join(base_path,relative_path)

