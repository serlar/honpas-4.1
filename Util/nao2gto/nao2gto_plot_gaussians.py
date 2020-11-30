#!/usr/bin/python3

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from math import exp
from ruamel.yaml import YAML

                    # ------------------------------------ #

# Calculate a Gaussian fit on an interval
def gaussfit(fit, rmin=0.0, rmax=5.0, points=50):

    xfit = np.linspace(rmin, rmax, points)
    yfit = np.zeros_like(xfit)

    for ir in range(len(xfit)):
        for ig in range(len(fit)):
            yfit[ir] += fit[ig][1] * exp(-fit[ig][0]*xfit[ir]*xfit[ir])

    return xfit, yfit

                    # ------------------------------------ #

# Select file to read
if ( len(sys.argv) >= 3 ):
    ref_path = sys.argv[1]
    fit_path = sys.argv[2]
else:
    fit_path = "nao2gto_ref.yml"
    fit_path = "nao2gto_fit.yml"

# Load NAO2GTO data
ref_yaml = YAML()
with open(ref_path, "r") as yml_file:
    ref_data = ref_yaml.load(yml_file)["orbitals"]
fit_yaml = YAML()
with open(fit_path, "r") as yml_file:
    fit_data = fit_yaml.load(yml_file)["orbitals"]

fit_lorb = [item["l"] for item in fit_data]
fit_zeta = [item["z"] for item in fit_data]

# Extract relevant data
nfns = len(fit_data)
cols = max(fit_lorb) - min(fit_lorb) + 1
rows = max(fit_zeta) - min(fit_zeta) + 1
print("Number of orbitals :", nfns)
print("Plot size          : {:d} rows, {:d} columns".format(rows,cols))
print("")

                    # ------------------------------------ #

# Build functions to plot
raw_vecs = []
fit_vecs = []
orb_sums = []
for ifns in range(nfns):

    # Take the input of GAUFRE as-is
    draw = fit_data[ifns]["raw"]
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

    # Calculate norms
    orb_sums.append(np.trapz(nraw, xraw))
    print("NORM(l={:d},z={:d}) = {:.7e}".format(
        fit_lorb[ifns], fit_zeta[ifns], orb_sums[-1]))
print("")

                    # ------------------------------------ #

# Plot orbitals
plt.rcParams["figure.figsize"] = [16, 4*rows+1]
fig, ax = plt.subplots(rows,cols)
plt.subplots_adjust(top=0.85, bottom=0.15, wspace=0.35, hspace=0.50)
fig.suptitle("Gaussian fitting of orbitals", fontsize=14)

for ifns in range(nfns):
    if ( fit_data[ifns]["do_fit"] ):
        fit_type = "normal"
    else:
        fit_type = "manual"

    icol = fit_lorb[ifns]
    irow = fit_zeta[ifns] - 1

    xraw, yraw = raw_vecs[ifns]
    xfit, yfit = fit_vecs[ifns]

    if ( rows == 1 ):
        bx = ax
    else:
        bx = ax[irow]

    bx[icol].plot(xraw, yraw, "bo", label="raw")
    bx[icol].plot(xfit, yfit, "r-", label="fit")

    bx[icol].yaxis.label.set_size(20)
    bx[icol].set(xlabel='r (Bohr)', ylabel='$\\frac{R_{nl}(r)}{r^l}$',
        title='L={:d}, $\zeta$={:d}, {:s} fit'.format(fit_lorb[ifns], fit_zeta[ifns], fit_type))
    bx[icol].grid()
    bx[icol].legend()

if ( (rows > 1) and (nfns % 2 == 1) ):
    ax[-1, -1].axis('off')

plt.show()
