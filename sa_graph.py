from __future__ import division
from numpy import *
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def sa_graph(start,sol):
    iterations, scores = loadtxt('start_' + str(start) + '_end_' + str(sol) + '.dat', unpack=True, usecols=[0,1])

    def movingaverage(values, window):
        weights = np.repeat(1.0, window) / window
        smas = np.convolve(values, weights, 'valid')
        return smas

    iterations = movingaverage(iterations, 10)
    scores = movingaverage(scores, 10)

    fig = plt.figure()

    ax = fig.add_subplot(111)
    ax.plot(iterations,scores,'k-',lw=1,label=str(start)+'_'+str(sol))
    ax.set_xlabel(r'Iteration')
    ax.set_ylabel(r'Score')
    plt.title('start_'+str(start) + '_end_'+str(sol))
    ax.figure.savefig('start_'+str(start)+'_end_'+str(sol)+'.pdf')
    plt.draw()





