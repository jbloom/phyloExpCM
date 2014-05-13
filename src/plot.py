"""Module for making plots with *matplotlib*.

Written by Jesse Bloom.


Dependencies
----------------
This module requires:

* *matplotlib*


Functions defined in this module
-----------------------------------

* *PlotLogLvsNParams* : plots log likelihood versus number of parameters.

* *PlotYearVersusDistance* : plots year separation versus distance

"""


import os
import math

# global variable _pylabavailable indicates if pylab / matplotlib present
try:
    import matplotlib
    matplotlib.use('pdf') # use PDF backend
    import pylab
    _pylabavailable = True
except ImportError:
    _pylabavailable = False


def PlotYearVersusDistance(plotfile, years, distances, xlabel='Year of isolation', ylabel='Distance', title=None):
    """Makes scatter plot of year versus distance.

    Requires ``matplotlib`` / ``pylab``; will raise an exception if
    not available.

    * *plotfile* : name of created PDF plot.

    * *years* : list of years, plotted on x-axis.

    * *distances* : list of distances, as same length as *years* and with
      corresponding entries. Plotted on y-axis.

    * *xlabel* : label for the x-axis, uses LaTex formatting.

    * *ylabel* : label for the y-axis, uses LaTex formatting.

    * *title* : if set to a non-empty value, is a string title placed
      above plot.
    """
    # check parameters
    if not _pylabavailable:
        raise ImportError("matplotlib and pylab are not available")
    if os.path.splitext(plotfile)[1].upper() != '.PDF':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    assert len(years) == len(distances), "years and distances not of same length"

    # plot size parameters
    (bigmargin, smallmargin) = (0.19, 0.02) # margins outside axes
    titlemargin = 0.15 # title height as fraction of plot
    (lmargin, rmargin, bmargin, tmargin) = (bigmargin, smallmargin, bigmargin, smallmargin)
    if title:
        tmargin += titlemargin
    plotmargin = 0.04 # add this much above / below last data point
    xsize = 2.1 # plot xsize in inches
    ysize = xsize * (1.0 - lmargin - rmargin) / (1.0 - tmargin - bmargin)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)

    # make plot
    figure = pylab.figure(figsize=(xsize, ysize), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.plot(years, distances, 'b.', markersize=5)
    (xmin, xmax, ymin, ymax) = (min(years), max(years), min(distances), max(distances))
    xmargin = plotmargin * (xmax - xmin)
    ymargin = plotmargin * (ymax - ymin)
    ax.set_xlim([xmin - xmargin, xmax + xmargin])
    ax.set_ylim([ymin - ymargin, ymax + ymargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    if title:
        pylab.title(title, size=10)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()




def PlotLogLvsNParams(plotfile, data_d):
    """Plots log likelihood (y-axis) versus number of parameters (x-axis).

    Requires *matplotlib* and *pylab*, and will raise an exception if 
    not available.

    CALLING VARIABLES:

    * *plotfile* is a string giving the name of the created plotfile. 
      Should end in the extension ``.pdf``.

    * *data_d* specifies the data points to be plotted. Each data point
      has three pieces of information:

        1) The *plotting_group* (a string).

        2) The number of parameters (number plotted on x-axis).

        3) The log likelihood (number plotted on y-axis).


      Points in the same *plotting_group* are plotted with the same symbol
      and have the same legend entry.

      *data_d* specifies these points as follows. *data_d* is a dictionary
      that is keyed by the *plotting_group* string. The entries for each
      string is a list of one or more 2-tuples. Each 2-tuple specifies
      a data point for that *plotting_group* as *(log_likelihood, nparams)*.

      Overall, *data_d* must specify at least two data points (although
      they can both be from the same plotting group if yo would like).

    """

    # check parameters
    if not _pylabavailable:
        raise ImportError("matplotlib and pylab are not available")
    if os.path.splitext(plotfile)[1].upper() != '.PDF':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if not (isinstance(data_d, dict) and data_d):
        raise ValueError("data_d is not a dictionary with at least one entries")

    # plot specifications
    markersize = 5 # size of symbols
    symbols = ['ro', 'bv', 'gs', 'y*', 'm+', 'k^', 'ch'] # list of plotting symbols 
    if len(data_d) > len(symbols):
        raise ValueError("Currently only enough symbols defined to for %d plotting groups" % (len(data_d)))
    (xsize, ysize) = (3.0, 3.0) # in inches
    (lmargin, rmargin, bmargin, tmargin) = (0.22, 0.08, 0.22, 0.08)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(xsize, ysize), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.xlabel('number of parameters', size=10)
    pylab.ylabel('log likelihood', size=10)
    plotlines = []
    igroup = 0
    symbols_names = []
    for (plotting_group, datalist) in data_d.iteritems():
        for (x, y) in datalist:
            linestyle = pylab.plot(x, y, symbols[igroup], markersize=markersize)
        symbols_names.append((plotting_group, linestyle[0]))
        igroup += 1
    pylab.legend([tup[1] for tup in symbols_names], [tup[0] for tup in symbols_names], loc='lower right', ncol=2, handlelength=1, columnspacing=0.75, labelspacing=0.8, handletextpad=0.25, numpoints=1)
    yformatter = matplotlib.ticker.ScalarFormatter()
    yformatter.set_scientific(True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
