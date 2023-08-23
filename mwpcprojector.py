# Short program to open and plot 2D histograms (in text form) alongside a circular window.
# It can calculate the window acceptance as a function of radius, outputting both a txt file and a graph.
#
#
# Jorge Romero 2020 joromero@jyu.fi

import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import font_manager
import os
import scienceplots

plt.style.use("science")
plt.rcParams["font.size"] = 50

np.seterr(all="ignore")


def averageSubtract(vector):
    avg = (vector.max() + vector.min()) / 2
    vector = vector - avg
    return vector


def Grain2DProjectionX(filename):
    x, _, n = np.genfromtxt(filename, unpack=True)

    resultx = []
    resultn = []
    localtotal = 0
    currentx = x[0]

    for i in range(len(x)):
        if currentx == x[i]:
            localtotal += n[i]
        else:
            resultx.append(currentx)
            resultn.append(localtotal)
            localtotal = n[i]
            currentx = x[i]

    resultx.append(currentx)
    resultn.append(localtotal)

    return resultx, resultn


def readGrain2DHistogram(filename):
    # Just read all three columns
    x, y, n = np.genfromtxt(filename, unpack=True)
    # print( x )
    # print( y )
    # print( n )

    # We assume that order of data is (x0,y0,n00), (x0,y1,n01), ...

    xlow = x[0]
    ylow = y[0]
    xhigh = x[-1]
    yhigh = y[-1]

    if usingDSSD:
        xavg = (xhigh + xlow) / 2
        yavg = (yhigh + ylow) / 2

        xlow -= xavg
        xhigh -= xavg
        ylow -= yavg
        yhigh -= yavg

    # This would fail if the order would be different
    ystep = y[1] - y[0]
    yn = int(np.round(float(yhigh - ylow) / ystep) + 1)

    # We know that total number of rows must be xn*yn. Solving xn.
    xn = int(n.size / yn)
    xstep = x[xn] - x[0]

    # Checking that dimensions are correct
    if xn * yn != n.size:
        raise Exception("Invalid file")

    # print( "xlow = {:}, xhigh = {:}, xstep = {:}, xn = {:}"
    #      .format( xlow, xhigh, xstep, xn ) )
    # print( "ylow = {:}, yhigh = {:}, ystep = {:}, yn = {:}"
    #     .format( ylow, yhigh, ystep, yn ) )

    # We have now information about the dimensions and we can
    # rearrange the shape of counts
    n = np.transpose(np.reshape(n, (xn, yn), order="C"))

    # print(n)

    return n, xlow, xhigh, ylow, yhigh


def readFileAndOutputArrays(filename):
    # reads the file and creates x-value, y-value and counts arrays. Also sums all of the counts in the file.

    xvalues, yvalues, nvalues = np.genfromtxt(filename, unpack=True)

    sum_all = np.sum(nvalues)
    return xvalues, yvalues, nvalues, sum_all


def plotSaveGrainHisto(name, xcent, ycent, WinRadius, Pdflag):
    filein = name + ".d2t"
    # read the file and output the correct count array and plot limits
    counts2d, xlow, xhigh, ylow, yhigh = readGrain2DHistogram(filein)

    # define the figure
    vertical = 2 if Pdflag else 1
    figsize = (30, 15.15) if Pdflag else (15.15, 15.3)

    fig, axes = plt.subplots(
        2,
        vertical,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    if Pdflag:
        mwpc1 = axes[0][0]
        proj1 = axes[1][0]
        mwpc2 = axes[0][1]
        proj2 = axes[1][1]
    else:
        mwpc1 = axes[0]
        proj1 = axes[1]

    # figure parameters
    mwpc1.set_ylabel("Y (mm)")
    # mwpc1.set_title(
    #    "$^{96}$Pd on MWPC at MARA Focal Plane. 64mm diameter window centred on ("
    #    + str(xcent)
    #    + ","
    #    + str(ycent)
    #    + ")."
    # )

    plt.rcParams["font.size"] = 50

    # fill the figure
    scale = 1.5 if Pdflag else 0.2
    mwpc1.imshow(
        counts2d, origin="lower", vmin=0, vmax=scale, extent=[xlow, xhigh, ylow, yhigh]
    )

    if Pdflag:
        mwpc2.imshow(
            counts2d,
            origin="lower",
            vmin=0,
            vmax=scale,
            extent=[xlow, xhigh, ylow, yhigh],
        )

    # draw two circles to represent the window sizes
    circ1 = Circle((xcent, ycent), radius=WinRadius, color=(1, 1, 1, 0.20))
    mwpc1.add_artist(circ1)
    mwpc1.set_ylim(ylow, yhigh + 0.25)

    if Pdflag:
        circ2 = Circle((-xcent + 2, ycent), radius=WinRadius, color=(1, 1, 1, 0.20))
        mwpc2.add_artist(circ2)
        mwpc2.set_ylim(ylow, yhigh + 0.25)

    # project on x axis
    if "Th" not in name:
        binwidth = 1.1
        countlabel = "Counts/0.25\,mm"
    else:
        binwidth = 4
        countlabel = "Counts/4\,mm"

    barcol = "dimgrey"
    projx, projn = Grain2DProjectionX(filein)
    proj1.bar(projx, projn, width=binwidth, color=barcol)
    if Pdflag:
        proj2.bar(projx, projn, width=binwidth, color=barcol)

    verticalcol = "C2"

    proj1.axvline(xcent - WinRadius, color=verticalcol, linewidth=5)
    proj1.axvline(xcent + WinRadius, color=verticalcol, linewidth=5)
    proj1.set_xlim(xlow, xhigh)
    proj1.set_ylim(0)

    if Pdflag:
        proj2.axvline(-xcent + 2 - WinRadius, color=verticalcol, linewidth=5)
        proj2.axvline(-xcent + 2 + WinRadius, color=verticalcol, linewidth=5)
        proj2.set_xlim(xlow, xhigh)
        proj2.set_ylim(0)
        mwpc2.set_yticks([])
        proj2.set_yticks([])
        proj2.set_xticks([-60, -30, 0, 30, 60])

    mwpc1.set_yticks([30, 15, 0, -15, -35])
    # proj1.set_yticks([])

    proj1.set_xticks([-60, -30, 0, 30, 60])

    padsize = 20 if Pdflag else 70
    proj1.set_ylabel(countlabel, labelpad=padsize)
    fig.supxlabel("X (mm)")

    # save the figure
    name = "./output" + name[1:] if Pdflag else "./output" + name[7:]
    fig.savefig(name + "_dual.pdf", transparent=True, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------

name = "./96Pd_Run20"  # no extension
window_centre_x = 9  # in mm
window_centre_y = 1  # in mm

"""
name = "./213Rn/213Rn_MWPC"  # no extension
window_centre_x = 3  # in mm
window_centre_y = 0  # in mm
"""

"""
name = "./226Th/226Th"  # no extension
window_centre_x = -24  # in mm
window_centre_y = 0  # in mm
"""

radius = 32  # in mm
max_radius = 40  # in mm
usingDSSD = False
Pdflag = True if "Pd" in name else False

plotSaveGrainHisto(name, window_centre_x, window_centre_y, radius, Pdflag)
