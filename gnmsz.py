#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

f = open("file.histo", "r")
kFrq = []
kCnt = []
for line in f :
    kFrq.append(int(line.split()[0]))
    kCnt.append(int(line.split()[1]))

# get local maxima in range 0 : 600
from scipy.signal import argrelextrema

maxima = argrelextrema(np.array(kCnt[20:1800]), np.greater)
mx = []
my = []
for maximum in maxima[0] :
    if kCnt[maximum+20] > 192000 :
        mx.append(maximum+21)
        my.append(kCnt[maximum+20])

globx = []
globy = []
for x, y in zip(mx, my) :
    print("Local maximum count: {} at {}X frequency".format(y, x))
    if len(globx) == 0 :
        globx.append(x) # Add the first value
        globy.append(y) # Add the first value
    elif globx[-1] != x and globx[-1]*1.1 < x : # if stored value of x + 10 percent < actual value of x --> we changed local peak
        globx.append(x)
        globy.append(y)
    elif globx[-1] != x and globx[-1]*1.1 >= x and globy[-1] < y : # if not > 10 % x value --> we still are in local peak and if current value > stored y value
        globx[-1] = x
        globy[-1] = y
print("Estimated max peaks at coverages: "+" ".join(str(x) for x in globx))

plt.style.use("seaborn")
fig, ax = plt.subplots(figsize=(15,10))

#ax.plot(kFrq[10:2000], kCnt[10:2000], marker=".", markerfacecolor="white")
ax.scatter(mx,my, marker="o", color="white", linewidth=0.5, edgecolor="blue", zorder=2)
ax.plot(kFrq[20:1000], kCnt[20:1000], zorder=1)
ax.set_xlabel("Coverage", fontsize=20)
ax.set_ylabel("$\it{k}$-mer count (frequency)", fontsize=20)

for x in globx :
    ax.plot([x,x], [0,1.05*max(kCnt[20:1000])], zorder=1, color="black", linewidth=0.7, linestyle="--")
    ax.annotate(xy=(x+5, 1.03*max(kCnt[20:1000])), s=str(x)+"X")

# CHANGE HERE DEPENDING ON PEAKS SHOWN
COV = 373

gnsz = sum([x*y for x,y in zip(kFrq[40:], kCnt[40:])])/COV

ax.annotate(xy=(-10, -0.04*max(kCnt[20:1000])), s="Estimated genome size: " + str(round(gnsz, 2))+"bp")
plt.show()

print("Estimated Genome Size: " + str(round(gnsz, 2)) + " (total sum divided by {})".format(COV))
