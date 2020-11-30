#!/usr/bin/python3

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from itertools import groupby
from ruamel.yaml import YAML

                    # ------------------------------------ #

# Calculate a Gaussian fit on an interval
def gaussfit(fit, rmin=0.0, rmax=5.0, points=50):

    xfit = np.linspace(rmin, rmax, points)
    yfit = np.zeros_like(xfit)

    for ir in range(len(xfit)):
        for ig in range(len(fit)):
            yfit[ir] += fit[ig][1] * np.exp(-fit[ig][0]*xfit[ir]*xfit[ir])

    return xfit, yfit

                    # ------------------------------------ #

# Select file to read
if ( len(sys.argv) >= 2 ):
    fit_path = sys.argv[1]
else:
    fit_path = "nao2gto_fit.yml"

# Load NAO2GTO data
fit_yaml = YAML()
with open(fit_path, "r") as yml_file:
    fit_data = fit_yaml.load(yml_file)["orbitals"]

fit_lorb = [item["qn_l"] for item in fit_data]
fit_zeta = [item["zeta"] for item in fit_data]
fit_orbs = [len(list(group)) for key, group in
    groupby([item["species"] for item in fit_data])]

# Extract relevant data
nfns = len(fit_data)
cols = max(fit_lorb) - min(fit_lorb) + 1
rows = (max(fit_zeta) - min(fit_zeta) + 1) * len(fit_orbs)
print("Number of orbitals :", nfns)
print("Plot size          : {:d} rows, {:d} columns".format(rows,cols))
print("")

                    # ------------------------------------ #

# Build functions to plot
raw_vecs = []
fit_vecs = []
orb_sums = []
fit_sigs = []
for ifns in range(nfns):

    # Take the input of GAUFRE as-is
    draw = fit_data[ifns]["points"]
    npts = len(draw)
    xraw = np.empty(npts)
    yraw = np.empty(npts)
    nraw = np.empty_like(xraw)

    # Transpose the existing data
    for ipts in range(npts):
        xraw[ipts] = draw[ipts][0]
        yraw[ipts] = draw[ipts][1]
        nraw[ipts] = (yraw[ipts]**2) * (xraw[ipts]**(2*(fit_lorb[ifns]+1)))
    raw_vecs.append([xraw, yraw])

    # Build a fitting curve
    fit_vecs.append(gaussfit(fit_data[ifns]["fit"]))
    sigmas = [1.0/np.sqrt(2*alpha) for alpha, beta in fit_data[ifns]["fit"]]
    fit_sigs.append(sigmas)

    # Calculate norms
    orb_sums.append(np.trapz(nraw, xraw))
    print("NORM(l={:d},z={:d}) = {:.7e}".format(
        fit_lorb[ifns], fit_zeta[ifns], orb_sums[-1]))
print("")

                    # ------------------------------------ #

# Plot orbitals
plt.rc('font', family='DejaVu Sans')
plt.rcParams["figure.figsize"] = [16, 4*rows+1]
fig, ax = plt.subplots(rows,cols)
plt.subplots_adjust(top=0.90, bottom=0.10, wspace=0.35, hspace=0.50)
fig.suptitle("Gaussian fitting of orbitals", fontsize=14)

for ifns in range(nfns):
    if ( fit_data[ifns]["do_fit"] ):
        fit_type = "yes"
    else:
        fit_type = "no"

    icol = ifns % cols
    irow = ifns // cols

    xraw, yraw = raw_vecs[ifns]
    xfit, yfit = fit_vecs[ifns]
    ymin = min(yraw)
    ymax = max(yraw)

    if ( rows == 1 ):
        bx = ax
    else:
        bx = ax[irow]

    bx[icol].plot(xraw, yraw, "bo", label="raw", markevery=10)
    bx[icol].plot(xfit, yfit, "r-", label="fit")
    bx[icol].axhline(color="black")
    bx[icol].vlines(fit_sigs[ifns], ymin, ymax, colors="green",
        linestyles="dashed")
    bx[icol].axvspan(max(xraw), 2.0*max(xraw), color="orange", alpha=0.7)

    bx[icol].set_xlim(left=0.0, right=1.05*max(xraw))
    bx[icol].set_ylim(bottom=(0.05-np.sign(ymin))*ymin,
        top=(0.05+np.sign(ymax))*ymax)
    bx[icol].yaxis.label.set_size(20)
    bx[icol].set(xlabel='r (Bohr)', ylabel='$\\frac{R_{nl}(r)}{r^l}$',
        title='{:s}: n={:d}, l={:d}, $\zeta$={:d}, fit: {:s} ({:d} GFs)'.format(
            fit_data[ifns]["species"], fit_data[ifns]["qn_n"],
            fit_lorb[ifns], fit_zeta[ifns], fit_type, len(fit_sigs[ifns])))
    bx[icol].grid()
    bx[icol].legend()

if ( (rows > 1) and (nfns % 2 == 1) ):
    ax[-1, -1].axis('off')

plt.show()
