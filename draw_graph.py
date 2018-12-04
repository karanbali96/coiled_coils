from __future__ import division
import numpy as np
from numpy import *
from scipy.interpolate import interp1d

def draw_graph(name):

    pH_x, charge_y = loadtxt(str(name)+'.dat',  unpack=True, usecols=[0,1])


    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    fig = plt.figure()

    ph_9 = np.interp(9, pH_x, charge_y)

    ax1 = fig.add_subplot(111)
    ax1.plot(pH_x, charge_y, 'k-', lw=1)
    ax2 = fig.add_subplot(111)
    ax2.axhline(y=ph_9, label='Charge at pH 9 is:'+str(ph_9), color='g')

    #f = interp1d(charge_y, pH_x)
    #pi = f(0)
    ph_5 = np.interp(5.8, pH_x, charge_y)

    ax3 = fig.add_subplot(111)
    ax3.axhline(y=ph_5, label= 'Charge at pH 5.8 is:'+str(ph_5), color='r')


    plt.legend(loc='best')
    plt.title(str(name))
    ax1.set_xlabel(r'pH')
    ax1.set_ylabel(r'Charge')
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax1.figure.savefig(str(name)+'.pdf')
    plt.draw()
