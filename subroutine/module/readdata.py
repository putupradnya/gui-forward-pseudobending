import numpy as np

def readEvent(path):
  data = open(path, 'r')
  file = data.readlines()

  for i in range(len(file)):
    file[i] = file[i].split()

  xevent = []
  yevent = []
  zevent = []

  for i in range(1, len(file)):
    xevent.append(np.float(file[i][0]))
    yevent.append(np.float(file[i][1]))
    zevent.append(np.float(file[i][2]))
  
  return xevent, yevent, zevent

def readStation(path):
  data = open(path, 'r')
  file = data.readlines()

  for i in range(len(file)):
    file[i] = file[i].split()

  xstation = []
  ystation = []
  zstation = []

  for i in range(1, len(file)):
    xstation.append(np.float(file[i][0]))
    ystation.append(np.float(file[i][1]))
    zstation.append(np.float(file[i][2]))
  
  return xstation, ystation, zstation

# xevent, yevent, zevent = readEvent('D:/Project/trainee/Main Program/hipocenter synthetic.dat')
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.scatter(xevent, yevent, zevent, marker='*', color='red', label='Source')
# plt.show()