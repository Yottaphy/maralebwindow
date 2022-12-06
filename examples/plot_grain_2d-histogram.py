import numpy as np
import scipy
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


def readGrain2DHistogram( filename ):
    # Just read all three columns
    x,y,n = scipy.genfromtxt( filename, unpack=True )
    print( x )
    print( y )
    print( n )

    # We assume that order of data is (x0,y0,n00), (x0,y1,n01), ...

    xlow  = x[0]
    ylow  = y[0]
    xhigh = x[-1]
    yhigh = y[-1]

    # This would fail if the order would be different
    ystep = y[1]-y[0]
    yn = int(np.round(float(yhigh-ylow)/ystep)+1)

    # We know that total number of rows must be xn*yn. Solving xn.
    xn = int(n.size / yn)
    xstep = x[xn]-x[0]

    # Checking that dimensions are correct
    if xn*yn != n.size:
        raise Exception("Invalid file")
    
    print( "xlow = {:}, xhigh = {:}, xstep = {:}, xn = {:}"
           .format( xlow, xhigh, xstep, xn ) )
    print( "ylow = {:}, yhigh = {:}, ystep = {:}, yn = {:}"
           .format( ylow, yhigh, ystep, yn ) )

    # We have now information about the dimensions and we can
    # rearrange the shape of counts
    n = np.transpose( np.reshape( n, (xn,yn), order='C' ) )

    print(n)

    return n, xlow, xhigh, ylow, yhigh



counts, xlow, xhigh, ylow, yhigh = readGrain2DHistogram( "test.d2t" )


fig = Figure( figsize=(15,10) )
canvas = FigureCanvas(fig)
ax     = fig.add_subplot(1,1,1)
ax.set_xlabel( 'X-axis title' )
ax.set_ylabel( 'Y-axis title' )
ax.set_title( 'Figure title' )

ax.imshow( counts, origin='lower', vmin=0, vmax=500,
           extent=[xlow,xhigh,ylow,yhigh])

fig.savefig('plot.png')
