
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  File: plots.py
#  Created: 04/26/18
#  Copyright 2018 Connor <mithrandir@hodor>
#
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

def plot_labels(type, ptype):
    if type == 1:
        title = 'Infinte Square Well'
    elif type == 2:
        title = 'Square Well'
    elif type == 3:
        title = 'Stepped Well'
    elif type ==4:
        title = 'Sloped Well'
    elif type == 5:
        title = 'Lattice'
    else:
        title = 'Simple Harmonic Osscilator'
    if ptype == 1:
        ylabel = "Psi"
    elif ptype == 2:
        ylabel = "Probablity Density (Psi^2)"
    return(title, ylabel)


def plot_setup(mtitle, ylabel, xMin, xMax):

    host = host_subplot(111, axes_class=AA.Axes)
    #plt.subplots_adjust(right=0.75)

    par1 = host.twinx()

    par1.axis["right"].toggle(all=True)

    host.set_xlim(xMin, xMax)
    plt.title(mtitle)
    host.set_xlabel("x - nm")
    host.set_ylabel(ylabel)
    par1.set_ylabel("Potential")

    return(host, par1)

def finish_plot(host, par1, p1, p2, y1min,y1max, y2min, y2max):
    host.set_ylim(y1min, y1max)
    par1.set_ylim(y2min, y2max)
    host.legend()

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())

    plt.draw()
    plt.show()
