#  -----------------------------------------------------------------------------------------------------
#
#  UKMON Orbit Plot
#
#  Plot UFO Orbit data for streams matching criteria.
#
#-- Acknowledgements
#
#   The R scripts in this module are an adaption of routines published originally as
#   Python scripts by Joe Hahn (http://gemelli,spacescience.org/~hahnjm)
#
#   Adapted by by Peter Campbell-Burns, 26 August 2016.
#
#  -----------------------------------------------------------------------------------------------------


Stream   = 'Perseids'

Auto_layout = False

x_limit = [-22,2]
y_limit = [-22,2]
z_limit = [-22,2]

#import modules used in code
import csv
import numpy as np
import matplotlib.pyplot as plt 
import math
import sys
from   kepler import *
import tkinter as tk
from   tkinter import filedialog
from   mpl_toolkits.mplot3d import Axes3D

#-- Function to plot 2D orbit

def orbit(a_element, e_element, i_element, o_element, w_element, Mp, col):
    Npoints = 10001
    GM = 1.0
    Mplot = np.linspace(0, 2*math.pi, num=Npoints)
    zero = 0.0
    max_error = 1.0e-10
    x_obj, y_obj, z_obj = el2xv(GM, a_element, e_element, zero, zero, zero, Mplot, max_error)
    axis.plot(x_obj, y_obj, z_obj, linestyle='-', color=col, linewidth=1)

#-- Function to plot 3D orbit (angles in Radians)

def orbit_3d(a_element, e_element, i_element, o_element, w_element, Mp, col):
    Npoints = 10001
    GM = 1.0
    zero = 0.0
    max_error = 1.0e-10
    Mplot = np.linspace(0, 2*math.pi, num=Npoints)
    x_obj, y_obj, z_obj = el2xv(GM, a_element, e_element, zero, zero, zero, Mplot, max_error)
    x_obj, y_obj, z_obj = Rz(x_obj, y_obj, z_obj, o_element)
    x_obj, y_obj, z_obj = Rx(x_obj, y_obj, z_obj, i_element)
    x_obj, y_obj, z_obj = Rz(x_obj, y_obj, z_obj, w_element)
    axis.plot(x_obj, y_obj, z_obj, linestyle='-', color=col, linewidth=0.5)       

#-- Define plot area

fig = plt.figure(figsize=(4, 4))
fig.suptitle('UFO Orbit plot')
axis = fig.add_subplot(111, projection='3d')
if Auto_layout : 
    plt.rcParams.update({'figure.autolayout': True})
else:
    axis.set_xlim(x_limit)
    axis.set_ylim(y_limit)
    axis.set_zlim(z_limit)
axis.azim = -60
axis.elev =  50
axis.set_xlabel('x (AU)')
axis.set_ylabel('y (AU)')
axis.set_zlabel('z (AU)')
plt.rcParams.update({'font.size': 9})
plt.rcParams.update({'figure.titlesize': 24})
axis.plot([0], [0], [0], marker='o', markersize=8.0, color='yellow')

#-- Plot Mercury orbit 

am     = 3.87098750E-01
em     = 2.05633707E-01
im_deg = 7.00427693E+00
om_deg = 4.83155706E+01
wm_deg = 2.91597880E+01
Mm_deg = 1.09448230E+02

am, em, im, om, wm, Mm = elements_deg2rad(am, em, im_deg, om_deg, wm_deg, Mm_deg)
orbit_3d(am, em, im, om, wm, Mm, 'darkgrey')

#-- Plot Venus Orbit

av = 0.72333199
ev = 0.00677323
iv_deg = 3.39471
ov_deg = 76.68069
wv_deg = 0.0
Mv_deg = 0.0

av, ev, iv, ov, wv, Mv = elements_deg2rad(av, ev, iv_deg, ov_deg, wv_deg, Mv_deg)
orbit_3d(av, ev, iv, ov, wv, Mv, 'darkgrey')

#-- Plot Earth Orbit

ae = 1.0
ee = 0.01671022
ie_deg = 0.00005
oe_deg = -11.26064
we_deg = 0.0
Me_deg = 0.0

ae, ee, ie, oe, we, Me = elements_deg2rad(ae, ee, ie_deg, oe_deg, we_deg, Me_deg)
orbit_3d(ae, ee, ie, oe, we, Me, 'green')

#-- File Dialog

root = tk.Tk()
root.withdraw()
UFO_FILE = filedialog.askopenfilename()

#-- Plot Meteor orbits

M_m_deg = 0.0
Line  = 0   
Selected = 0

with open(UFO_FILE, 'rt') as f:
    reader = csv.reader(f) 

    # locate data in CSV
    h = next(reader)
    l = len(h)
    for i in range(0,l):
        if h[i] == "_a":
            _a = i
        elif h[i] == "_q":
            _q = i
        elif h[i] == "_e":
            _e = i
        elif h[i] == "_p":
            _p = i           
        elif h[i] == "_peri":
            _w = i
        elif h[i] == "_node":
            _o = i
        elif h[i] == "_incl":
            _i = i          
    if not ('_a' in locals() and '_q' in locals() and '_e' in locals() and '_p' in locals() and '_w' in locals() and '_o' in locals() and '_e' in locals()): 
        print('ERROR - Invalid CSV')        
        sys.exit()
        
    # Create plot for each row    
    for row in reader:
        if row[1].strip() == Stream:
            Col = 'blue'
            Selected = Selected + 1
            
            # Extract UFO dta
            a_m     = float(row[_a]) # Semi mahor axix
            q_m     = float(row[_q]) 
            e_m     = float(row[_e]) # Eccentricity
            p_m     = float(row[_p])
            w_m_deg = float(row[_w])
            o_m_deg = float(row[_o])
            i_m_deg = float(row[_i])
            
            a_m, e_m, i_m, o_m, w_m, M_m = elements_deg2rad(a_m, e_m, i_m_deg, o_m_deg, w_m_deg, M_m_deg)
            orbit_3d(a_m, e_m, i_m, o_m, w_m, M_m, Col)
        Line = Line + 1

print("Number of meteors plotted: ",Selected)
 
axis.text3D(-0.13, 0.54, -0.05, 'M', color='darkgrey')
axis.text3D(0.0, 0.88, -0.05, 'V', color='darkgray')
axis.text3D(0.15, 1.15, -0.05, 'E', color='green')
axis.text3D(0.7, -3.0, 0.5, 'Geminids', color='red')

plt.show(block=False)
plt.draw()




