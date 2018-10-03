#!/usr/bin/python
import subprocess
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import re
import cPickle
import numpy
from math import log
from math import sqrt


FIRST_PHINDEX = 3
LAST_PHINDEX = 23
window_size = 65536
eps = 2**(-8)

def plot_memory(memory):
    average_memory= memory
    deltas = [0.01, 0.05, 0.1, 0.2] #actually its counters parameters (counters = 1/epsilon)
    MS = 12
    LW = 4
    print deltas
    print average_memory["ECM"]
    plt.plot(deltas, average_memory["ECM"],"-D",label="$ECM$", markersize=MS, linewidth=LW, c="darkcyan")

    #plt.xscale("log",basex=2)
    #plt.yscale("log",basey=2)

    ##for algorithm in speed:
    ##    plt.plot(phis, speed[algorithm],label=algorithm)
    #plt.xscale("log",basex=2)
    #plt.yscale("log",basex=2)

    #ticks,labels = plt.xticks()
    #plt.xticks(ticks[::2],labels[::2])
    plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))

    plt.xlabel("Error Probability $\delta$", fontsize=36)
    ylabel_str = "Memory [MB]"
    plt.ylabel(ylabel_str, fontsize=36)
    plt.tick_params(labelsize=20)
    plt.xlim(0, 0.2)
    plt.ylim(16, 30)
    plt.xticks([0.01, 0.05, 0.1, 0.2], ('0.01%', '0.05%', '0.1%', '0.2%'))

    plt.legend(loc="best") # keys of the graphs
    plt.tight_layout()
    plt.savefig('test_delta_memory.png')
    plt.clf()

def calc_memory(delta):
    RSS_WCSS = 108 * (4/eps)
    RSS_HIT = 108 * (5/eps)
    RSS_ACC1 = 108 * (5/eps)

    B_WCSS = (4/eps) * (4 + ((log((4/eps), 2))/8))

    b_WCSS = 4 * (4/eps)

    ACC1 = (5/eps) * (5/eps) * (4 + log(5/eps, 2)/8)
    ACC2 = (5/eps) * 2 * sqrt((5/eps)) * (4 + log(5/eps, 2)/8)
    acc = []
    for i in range(1, 32):
        acc.append((5/eps)**(1 + (1.0/ i)) * i  * (4 + log(5/eps, 2)/8) + RSS_ACC1)
    HIT = (5/eps) * log(5/eps, 2) * (4 + log(5/eps, 2)/8)

    WCSS = RSS_WCSS + B_WCSS + b_WCSS
    ACC1 = RSS_ACC1 + ACC1
    ACC2 = RSS_ACC1 + ACC2
    HIT = RSS_ACC1  + HIT
    RAW = (4/eps) * WCSS
    oneOverDelta = 1/(delta/100);
    logOneOverDelta = log(oneOverDelta,2);
    memory["ECM"].append((logOneOverDelta * eps**(-2) * 32)/ 1e6)

memory = dict()
memory["ECM"]=[]
deltas = [0.01, 0.05, 0.1, 0.2] #actually its counters parameters (counters = 1/epsilon)

for delta in deltas:
    out = calc_memory(delta)
plot_memory(memory)
