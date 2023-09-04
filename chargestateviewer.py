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

# set font

plt.style.use("science")
plt.rcParams["font.size"] = 40


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


def plotSaveGrainHisto(name):
    filein = name + ".d2t"
    # read the file and output the correct count array and plot limits
    counts2d, xlow, xhigh, ylow, yhigh = readGrain2DHistogram(filein)

    # define the figure
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 16.16),
        sharex=False,
        gridspec_kw={"hspace": 0, "wspace": 0},
    )
    mwpc = axes[1]
    colb = axes[0]
    colb.set_visible(False)

    # figure parameters
    mwpc.set_ylabel("Y (mm)")
    mwpc.set_xlabel("X (mm)")
    # mwpc.set_title(
    #    "$^{96}$Pd on MWPC at MARA Focal Plane. 64mm diameter window centred on ("
    #    + str(xcent)
    #    + ","
    #    + str(ycent)
    #    + ")."
    # )

    # fill the figure
    image = mwpc.imshow(
        counts2d,
        cmap="Greys",
        origin="lower",
        vmin=0,
        vmax=5000,
        extent=[xlow, xhigh, ylow, yhigh],
    )

    # Text
    mwpc.text(-65, -30, "q=29$^+$", size="small")
    mwpc.text(-40, -30, "q=28$^+$", size="small")
    mwpc.text(-15, -30, "q=27$^+$", size="small")
    mwpc.text(10, -30, "q=26$^+$", size="small")
    mwpc.text(35, -30, "q=25$^+$", size="small")

    # colourbar
    clb = fig.colorbar(image, ax=colb, orientation="horizontal", aspect=50)
    clb.set_label("Counts")
    clb.ax.xaxis.set_ticks_position("top")
    clb.ax.xaxis.set_label_position("top")

    # save the figure
    fig.savefig(name + "_MWPCimage.pdf", transparent=True, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------


name = "./M17_MWPC"  # no extension

usingDSSD = False

plotSaveGrainHisto(name)
