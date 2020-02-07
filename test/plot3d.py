from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# create a 21 x 21 vertex mesh
model = np.load('D:/Project/trainee/ModelVelocity.npy')

cbar = np.reshape(model, 40*40*40)
x =  np.zeros(40) + 1000 #np.r_[0:20000:500]
y = np.r_[0:20000:500]
z = np.r_[0:20000:500]

plt.close('all')
fig = plt.figure()
ax = fig.gca(projection='3d')
mappable = plt.cm.ScalarMappable(cmap=plt.cm.jet)
mappable.set_array(model)

X, Y = np.meshgrid(x, y)

# Vertical velocity
data = np.zeros((40,40))
for i in range(40):
    dat = model[:][i][0]
    for j in range(40):
        data[j][i] = dat[j]

Z = data #model[:][:][0]
Z1 = model[:][:][7]
Z2 = model[:][:][19]
# cset = [[],[],[]]

# this is the example that worked for you:
# ax.contourf(X, Y, Z, zdir='z',cmap='jet')

# now, for the x-constant face, assign the contour to the x-plot-variable:
ax.contourf(Z, Y, X, zdir='x',cmap='jet')

# likewise, for the y-constant face, assign the contour to the y-plot-variable:
# ax.contourf(X, Z2, Y, zdir='y',cmap='jet')
ax.invert_zaxis()
ax.set_xlabel('X Coor')
ax.set_ylabel('Y Coor')
ax.set_zlabel('Z Coor')

plt.colorbar(mappable)
plt.show()
