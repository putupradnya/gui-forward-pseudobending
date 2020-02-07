import numpy
from mayavi import mlab
import numpy as np
# from enthought.mayavi.modules.scalar_cut_plane import ScalarCutPlane

# create a 21 x 21 vertex mesh
model = np.load('D:/Project/trainee/ModelVelocity.npy')

# Vertical velocity
data = np.zeros((40,40))
for i in range(40):
    dat = model[:][i][0]
    for j in range(40):
        data[j][i] = dat[j]


# creating volume that increases in value
img3d = np.arange(20)
img3d = np.expand_dims(img3d, axis=1)
img3d = np.expand_dims(img3d, axis=2)
img3d = np.tile(img3d, (1, 20, 20))

fig = mlab.figure()
src = mlab.pipeline.scalar_field(img3d)

# Plotting two cut planes
cp2 = mlab.pipeline.scalar_cut_plane(src, plane_orientation='y_axes')
cp2.implicit_plane.widget.enabled = False
cp3 = mlab.pipeline.scalar_cut_plane(src, plane_orientation='z_axes')
cp3.implicit_plane.widget.enabled = False
mlab.view(azimuth=50, elevation=None)
mlab.outline()
mlab.scalarbar(orientation='vertical')
mlab.show()
# obj = mlab.imshow(data)
# mlab.show()