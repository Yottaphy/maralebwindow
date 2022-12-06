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

# set font

plt.rcParams["font.family"] = "Latin Modern Roman"
plt.rcParams["font.size"] = 50


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


def countInWindow(xvalues, yvalues, nvalues, sum_all, xcent, ycent, radius):
    # sums all of the counts within the circular window defined by a specified radius and center point (xcent, ycent) and calculates the ratio of those counts to the total counts

    xrec = list(map(lambda x: (x - xcent) ** 2, xvalues))
    yrec = list(map(lambda x: (x - ycent) ** 2, yvalues))

    mask = np.sqrt(np.add(xrec, yrec)) < radius

    sum_selected = np.sum(nvalues[mask])

    ratio = sum_selected / sum_all

    return sum_selected, ratio


def multiRadiusPlot(name, xcent, ycent, max_radius, selectedRadius):
    # This function plots and writes a txt file for multiple radii
    filename = name + ".d2t"
    # reads the filename, feeding arrays of x values, y values, number of counts and the total number of counts to the window counter
    xvalues, yvalues, nvalues, sum_all = readFileAndOutputArrays(filename)

    if usingDSSD:
        xvalues = averageSubtract(xvalues)
        yvalues = averageSubtract(yvalues)

    r = []
    s = []
    # open file for writing
    f = open(name + "_YieldvRadius_" + str(xcent) + "_" + str(ycent) + ".txt", "w+")
    f.write(
        "Window centred at x = " + str(xcent) + " mm and y = " + str(ycent) + " mm\n"
    )
    f.write("\n Radius (mm)\tAcceptance (%)\n")

    # calculates the sum of counts within windows of variable width and the ratio of that to total counts in the histogram
    for i in range(0, max_radius + 1):
        sum_s, ratio = countInWindow(
            xvalues, yvalues, nvalues, sum_all, xcent, ycent, i
        )
        r.append(ratio * 100)
        if i == selectedRadius:
            print(
                "Transmission:",
                ratio * 100,
                "% for a",
                selectedRadius,
                "mm radius window.",
            )
        s.append(sum_s)
        f.write(
            str(i) + "\t" + str(round(ratio * 100, 3)) + "\n"
        )  # writes the radius in mm and the percentage acceptance

    f.close()  # closes the file after the writing

    # plots ratio vs window radius
    plt.plot(r)
    plt.ylabel("Counts within window (%)")
    plt.xlabel("Window Radius (mm)")
    # plt.title(
    #    "Percentage acceptance of the gas cell window\nfor different radii when centering on ("
    #    + str(xcent)
    #    + ","
    #    + str(ycent)
    #    + ")."
    # )
    plt.grid(which="major", axis="both")
    plt.savefig(name + "_multiradius.pdf")


def plotSaveGrainHisto(name, xcent, ycent, WinRadius):
    filein = name + ".d2t"
    # read the file and output the correct count array and plot limits
    counts2d, xlow, xhigh, ylow, yhigh = readGrain2DHistogram(filein)

    # define the figure
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 16.16),
        sharex=True,
        gridspec_kw={"hspace": 0, "wspace": 0},
    )
    mwpc = axes[0]
    proj = axes[1]

    # figure parameters
    mwpc.set_ylabel("Y (mm)")
    # mwpc.set_title(
    #    "$^{96}$Pd on MWPC at MARA Focal Plane. 64mm diameter window centred on ("
    #    + str(xcent)
    #    + ","
    #    + str(ycent)
    #    + ")."
    # )

    # fill the figure
    scale = 1.5 if "Pd" in name else 0.2
    mwpc.imshow(
        counts2d, origin="lower", vmin=0, vmax=scale, extent=[xlow, xhigh, ylow, yhigh]
    )

    # draw two circles to represent the window sizes
    circ1 = Circle((xcent, ycent), radius=WinRadius, color=(1, 1, 1, 0.20))
    mwpc.add_artist(circ1)
    mwpc.set_ylim(ylow, yhigh + 0.25)

    # project on x axis
    binwidth = 1.1 if "Th" not in name else 4

    projx, projn = Grain2DProjectionX(filein)
    proj.bar(projx, projn, width=binwidth, color="dimgrey")

    proj.axvline(xcent - WinRadius, color="C1", linewidth=5)
    proj.axvline(xcent + WinRadius, color="C1", linewidth=5)

    proj.set_xlim(xlow, xhigh)
    proj.set_ylim(0, 5.5)
    proj.set_ylabel("Counts")
    proj.set_xlabel("X (mm)")

    # save the figure
    fig.savefig(name + "_map.pdf", transparent=True, bbox_inches="tight")
    fig.savefig(name + "_map.png", bbox_inches="tight")


# ----------------------------------------------------


"""
name = "./96Pd_Run20"  # no extension
window_centre_x = 9  # in mm
window_centre_y = 1  # in mm
"""

"""
name = "./213Rn/213Rn_MWPC"  # no extension
window_centre_x = 3  # in mm
window_centre_y = 0  # in mm

"""
name = "./226Th/226Th"  # no extension
window_centre_x = -24  # in mm
window_centre_y = 0  # in mm

radius = 32  # in mm
max_radius = 40  # in mm
usingDSSD = False

plotSaveGrainHisto(name, window_centre_x, window_centre_y, radius)
# multiRadiusPlot(name, window_centre_x, window_centre_y, max_radius, radius)
