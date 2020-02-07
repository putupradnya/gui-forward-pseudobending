import numpy as np
import subroutine.module.basic as bs

def forwardProcess(self,niter,cacah,split,x,y,z,vxyz,xs,ys,zs,xr,yr):
  dvx, dvy, dvz = bs.diff_cube(vxyz, x, y, z)

  xb, yb, zb = bs.iterbending(niter,cacah,x,y,z,vxyz,xs,ys,zs,xr,yr,dvx,dvy,dvz)

  time = []

  # Calculate ray propagation time
  for i in range(len(xb) - 1):
    point_x = np.linspace(xb[i], xb[i+1], split)
    point_y = np.linspace(yb[i], yb[i+1], split)
    point_z = np.linspace(zb[i], zb[i+1], split)

    for k in range(len(point_x) - 1):
      length_seg = bs.distance(point_x[k], point_y[k], point_z[k], point_x[k+1], point_y[k+1], point_z[k+1])
      xmean = (point_x[k] + point_x[k+1])/2
      ymean = (point_y[k] + point_y[k+1])/2
      zmean = (point_z[k] + point_z[k+1])/2
      a, b, c = bs.index(x, y, z, xmean, ymean, zmean)
      time.append(length_seg/vxyz[a][b][c])

  time_total = np.sum(time)
  # print('time elapsed: ', time_total)

  return xb, yb, zb, time_total