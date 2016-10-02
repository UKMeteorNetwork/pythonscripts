#kepler.py
#
#  by Joe Hahn, jhahn@spacescience.org, 26 January 2014.
#  Adapted by by Peter Campbell-Burns,  26 August 2016.

import numpy as np
import math

def keplers_eqn(e, M, max_error):
    #Solve kepler's equation via Danby's algorithm.
    twopi = 2*math.pi
    M_mod = np.fmod(M, twopi)
    arr = np.array([1.0, 2.0])
    if (type(M_mod) != type(arr)): M_mod=np.array([M_mod])
    sgn = M_mod*0 + 1.0
    sgn[np.sin(M_mod) < 0.0] = -1.0
    E = M_mod + 0.85*sgn*e
    max_iteration = 15
    for i in range(max_iteration):   
        es = e*np.sin(E)
        ec = e*np.cos(E)
        f = E - es - M_mod
        error = np.fabs(f).max()
        if error < max_error: break
        df = 1.0 - ec
        ddf = es
        dddf = ec
        d1 = -f/df
        d2 = -f/(df + d1*ddf/2.0)
        d3 = -f/(df + d2*ddf/2.0 + d2*d2*dddf/6.0)
        E = E + d3
    if error > max_error: 
        print("***Warning*** keplers_eqn() failed to converge, error = ", error)
    return E


def el2polar(GM, a, e, M, max_error):
    #convert orbit elements a,e,M to polar coordinates r,f,z,vr,vt,vz in the z=0 orbit plane.
    E = keplers_eqn(e, M, max_error)
    r = a*(1 - e*np.cos(E))
    f = 2.0*np.arctan( np.sqrt((1 + e)/(1 - e))*np.tan(E/2) )
    n = np.sqrt(GM/a/a/a)
    vr = e*a*n*np.sin(f)/np.sqrt(1 - e*e)
    vt = a*n*(1 + e*np.cos(f))/np.sqrt(1 - e*e)
    z = vz = r*0.0
    return r, f, z, vr, vt, vz

def polar2cartesian(r, theta, vr, vt):
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    vx = vr*np.cos(theta) - vt*np.sin(theta)
    vy = vr*np.sin(theta) + vt*np.cos(theta)
    return x, y, vx, vy

def el2xv(GM, a, e, I, O, w, M, max_error):
    #convert orbit elements to cartesian coordinates and velocities.
    r, f, z, vr, vt, vz = el2polar(GM, a, e, M, max_error)
    x, y, vx, vy = polar2cartesian(r, f, vr, vt)
    x, y, z = Rz(x, y, z, -w)
    vx, vy, vz = Rz(vx, vy, vz, -w)
    x, y, z = Rx(x, y, z, -I)
    vx, vy, vz = Rx(vx, vy, vz, -I)
    x, y, z = Rz(x, y, z, -O)
    vx, vy, vz = Rz(vx, vy, vz, -O)
    return x, y, z, vx, vy, vz

def xv2el(GM, x, y, z, vx, vy, vz):
    #convert cartesian coordinates and velocities to orbit elements.
    hx = y*vz - z*vy
    hy = z*vx - x*vz
    hz = x*vy - y*vx
    hxy2 = hx*hx + hy*hy
    hxy = np.sqrt(hxy2)
    h2 = hxy2 + hz*hz
    h = np.sqrt(h2)
    I = np.arccos(hz/h)
    r = np.sqrt(x*x + y*y + z*z)
    ex = (vy*hz - vz*hy)/GM - x/r
    ey = (vz*hx - vx*hz)/GM - y/r
    ez = (vx*hy - vy*hx)/GM - z/r
    e = np.sqrt(ex*ex + ey*ey + ez*ez)
    O = np.arctan2(hx, -hy)
    ec = ex*np.cos(O) + ey*np.sin(O)
    es = ez/np.sin(I)
    v2 = vx*vx + vy*vy + vz*vz
    a = 1.0/(2.0/r - v2/GM)
    w = np.arctan2(es, ec)
    rv = x*vx + y*vy + z*vz
    ec = 1.0 - r/a
    es = rv/np.sqrt(GM*a)
    E = np.arctan2(es, ec)
    M = E - es
    return a, e, I, O, w, M

def Rx(x, y, z, angle):
    #rotate coordinate system about the x axis by angle, in radians.
    s = np.sin(angle)
    c = np.cos(angle)
    xr = x
    yr = y*c + z*s
    zr =-y*s + z*c
    return xr, yr, zr

def Ry(x, y, z, angle):
    #rotate coordinate system about the y axis by angle, in radians.
    s = np.sin(angle)
    c = np.cos(angle)
    xr = x*c - z*s
    yr = y
    zr = x*s + z*c
    return xr, yr, zr

def Rz(x, y, z, angle):
    #rotate coordinate system about the z axis by angle, in radians.
    s = np.sin(angle)
    c = np.cos(angle)
    xr = x*c + y*s
    yr =-x*s + y*c
    zr = z
    return xr, yr, zr

def elements_deg2rad(a, e, I, O, w, M):
    #convert the angular orbit elements from degrees to radians.
    f = math.pi/180
    return a, e, I*f, O*f, w*f, M*f


