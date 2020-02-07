import numpy as np

def diff(x, y):
    yb2 = []
    for i in range(len(y)):
        if (i <= (len(y) - 2)):
            yb2.append((y[i + 1] - y[i]) / ((x[i + 1] - x[i])))
        else:
            yb2.append((y[i] - y[i - 1]) / ((x[i] - x[i - 1])))
    return yb2

def find(z, a, dvzi):
    for i in range(len(z)):
        if (z[i] <= a):
            x = dvzi[i]
        else:
            break
    return x


def mid(x1, x2):
    xmid = 0.5 * (x1 + x2)
    return xmid


def diff_cube(vxyz, x, y, z):
    dvx = np.zeros(vxyz.shape)
    for i in range(vxyz.shape[2]):
        for j in range(vxyz.shape[1]):
            dvx[:, j, i] = diff(x, vxyz[j, i, :])
    dvy = np.zeros(vxyz.shape)
    for i in range(vxyz.shape[2]):
        for j in range(vxyz.shape[0]):
            dvy[j, :, i] = diff(y, vxyz[i, :, j])
    dvz = np.zeros(vxyz.shape)
    for i in range(vxyz.shape[0]):
        for j in range(vxyz.shape[1]):
            dvz[i, j, :] = diff(z, vxyz[:, j, i])
    return dvx, dvy, dvz


def index(x, y, z, xi, yi, zi):
    i = np.abs(np.fix((xi - x[0]) / (x[1] - x[0])))
    j = np.abs(np.fix((yi - y[0]) / (y[1] - y[0])))
    k = np.abs(np.fix((zi - z[0]) / (z[1] - z[0])))
    return int(i), int(j), int(k)

def ray_bending(x1, y1, z1, x2, y2, z2, dvx, dvy, dvz, xmid, ymid, zmid, dvxm, dvym, dvzm, Vmid, V1, V2):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    L2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
    nL = dx * dvx + dy * dvy + dz * dvz
    nx = dvx - nL * dx / L2
    ny = dvx - nL * dy / L2
    nz = dvz - nL * dz / L2

    nl = np.sqrt(pow(nx, 2) + pow(ny, 2) + pow(nz, 2))
    nx = nx / (nl+10**(-6))
    ny = ny / (nl+10**(-6))
    nz = nz / (nl+10**(-6))

    l = pow(x2 - xmid, 2) + pow(y2 - ymid, 2) + pow(z2 - zmid, 2)
    c = 0.5 / V1 + 0.5 / V2
    ra = pow(0.25 * (c * Vmid + 1) / (c * (nx * dvxm + ny * dvym + nz * dvzm)+10**(-6)), 2)
    rb = 0.5 * l / (c * Vmid)
    Rc = -0.25 * (c * Vmid + 1) / (c * (nx * dvxm + ny * dvym + nz * dvzm)+10**(-6)) + np.sqrt(ra + rb)
    xnew = xmid + nx * Rc
    ynew = ymid + ny * Rc
    znew = zmid + nz * Rc
    return xnew, ynew, znew

def iterbending(niter,cacah,x,y,z,vxyz,xs,ys,zs,xr,yr,dvx,dvy,dvz):
    xi = [xr, xs]
    yi = [yr, ys]
    zi = [0, zs]

    xi = np.linspace(xi[0], xi[1], cacah + 1)
    yi = np.linspace(yi[0], yi[1], cacah + 1)
    zi = np.linspace(zi[0], zi[1], cacah + 1)

    xb = np.zeros(len(xi))
    yb = np.zeros(len(xi))
    zb = np.zeros(len(xi))

    for j in range(niter):
        for i in range(cacah - 1):
            x1 = xi[i]
            y1 = yi[i]
            z1 = zi[i]
            x2 = xi[i + 2]
            y2 = yi[i + 2]
            z2 = zi[i + 2]
            xmid = mid(x1, x2)
            ymid = mid(y1, y2)
            zmid = mid(z1, z2)
            o, p, q = index(x, y, z, x1, y1, z1)
            if o < len(dvx[:,0,0]) and p < len(dvx[0,:,0]) and q < len(dvx[0,0,:]):
                dvxi = dvx[o, p, q]
                dvyi = dvy[o, p, q]
                dvzi = dvz[o, p, q]
                V1 = vxyz[o, p, q]
            else:
                return xi, yi, zi
            o, p, q = index(x, y, z, x2, y2, z2)
            if o < len(dvx[:,0,0]) and p < len(dvx[0,:,0]) and q < len(dvx[0,0,:]):
                V2 = vxyz[o, p, q]
            else:
                return xi, yi, zi
            o, p, q = index(x, y, z, xmid, ymid, zmid)
            if o < len(dvx[:,0,0]) and p < len(dvx[0,:,0]) and q < len(dvx[0,0,:]):
                dvxm = dvx[o, p, q]
                dvym = dvy[o, p, q]
                dvzm = dvz[o, p, q]
                Vmid = vxyz[o, p, q]
            else:
                return xi, yi, zi
            xn, yn, zn = ray_bending(x1, y1, z1, x2, y2, z2, dvxi, dvyi, dvzi, xmid, ymid, zmid, dvxm, dvym, dvzm,
                                         Vmid,V1, V2)
            if np.isnan(xn) == False and np.isnan(yn) == False and np.isnan(zn) == False:
                o, p, q = index(x, y, z, xn, yn, zn)
                if o < len(dvx[:,0,0]) and p < len(dvx[0,:,0]) and q < len(dvx[0,0,:]):
                    xb[i + 1] = xn
                    yb[i + 1] = yn
                    zb[i + 1] = zn
                else:
                    return xi, yi, zi
            else:
                return xi, yi, zi

        xb[0] = xi[0]
        yb[0] = yi[0]
        zb[0] = zi[0]
        xb[-1] = xi[-1]
        yb[-1] = yi[-1]
        zb[-1] = zi[-1]
        if min(zb) >= min(zi) and max(zb) <= max(zi):
            xi = xb
            yi = yb
            zi = zb
        else:
            return xi, yi, zi
    return xi,yi,zi

def distance(x1, y1, z1, x2, y2, z2):
  length = np.sqrt(pow((x2- x1),2) + pow((y2 - y1), 2) + pow((z2-z1),2))
  return length