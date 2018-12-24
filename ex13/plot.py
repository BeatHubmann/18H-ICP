#!/usr/bin/env python3

# based on and adapted from helper code provided by HPCSEI course team in HS18 @ETHZ

import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

import glob
import os

datafiles = sorted(glob.glob('wave_*.dat'))
print("Plotting {} frame(s).".format(len(datafiles)))
print("This will take about {}MB of space!".format(int(1.23 * len(datafiles))))
for filename in datafiles:
    print(filename)
    L = np.loadtxt(filename, skiprows=1)
    x = L[:, 0]
    y = L[:, 1]

    plt.clf()
    ax = plt.gca()
    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(-1.0, 1.0)
    plt.plot(x, y)
    # vmin and vmax is the range of phi values.
    plt.savefig(filename + '.raw')

print("Running ffmpeg to generate a video...")
os.system('ffmpeg -v panic -stats -framerate 20 -pix_fmt rgba -s 640x480 -i "wave_%05d.dat.raw" -c:v libx264 -r 20 -pix_fmt yuv420p -y "wave.mp4"')
print("Done!")
