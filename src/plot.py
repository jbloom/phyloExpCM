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

* *PlotSiteLikelihoods* : plots site likelihood comparison

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


def PlotSiteLikelihoods(sites, classifications, classificationtypes, xvalues, xlabel, ylabel, plotfile, symmetrize_axis=True, title='', alpha=1):
    """Plots comparison of site likelihoods.

    Requires ``matplotlib`` / ``pylab``, will raise an exception if
    not available.

    This function plots a comparison of the log likelihoods at a set of sites.

    Typically you will have computed the difference in the log likelihoods for each
    site over two different models. You may then have classified the sites into
    one or more categories. This function plots points showing the differences
    in log likelihoods for all sites in each classification.

    The each set of data is plotted using a box plots.

    All strings (*xlabel*, *ylabel*, *classificationtypes*, *title*) use ``LaTex``
    formatting.

    CALLING VARIABLES:

    * *sites* is a list of integer sites (should all be unique)

    * *classifications* is a dictionary that is keyed by all integers in *sites*.
      The values should be some set of (relatively few) strings that define the site
      classification. Don't make too many classifications, as a separate row is plotted
      for each.

    * *classificationtypes* is a list of the classifications (the unique values in
      *classifications*) in the order they should be plotted from bottom to top on the
      y-axis.

    * *xvalues* is a dictionary keyed by all integers in *sites*, with the value being
      the number plotted on the x-axis. Typically this would be the difference in
      log likelihoods between model 1 and model 2.

    * *xlabel* is the string label for the x-axis. This might be something like
      *LogL for model 1 - LogL for model 2*.

    * *ylabel* is the string label for the y-axis.

    * *plotfile* is the name of the created plot. Must have the extension
      ``.pdf``.

    * *symmetrize_axis* is a Boolean switch specifying if we make the x-axis
      range equal at both ends (i.e. goes from -10 to 10). Is *True* by default.

    * *title* is an optional string above the plot. Is empty (*''*) by default.
      Uses ``LaTex`` formatting.

    * *alpha* is the transparency of the plotted points. An *alpha* of 1 means
      no transparency, and *alpha* less than one means transparent. Is 1 by
      default.
    """
    if not _pylabavailable:
        raise ImportError("This function requires matplotlib / pylab")
    if os.path.splitext(plotfile)[1].upper() != '.PDF':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    # set up plot
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    xdatamargin = 1.08 # make data margins this many times greater than maximal point
    xsize = 2.0 # plot width in inches
    lmargin = 0.95 # left margin width in inches
    rmargin = 0.05 # right margin width in inches
    bmargin = 0.6 # bottom margin width in inches
    tmargin = 0.05 # top margin width in inches
    plotheight = 0.42 # height of plot used for each classification
    if title:
        titlemargin = 0.25 # title margin width in inches
    else:
        titlemargin = 0.0 # title margin width in inches
    classificationvalues = dict([(classification, True) for classification in classifications.itervalues()]).keys()
    assert len(classificationvalues) == len(classificationtypes), "Not the same number of values in classifications and entries in classificationtypes"
    for c in classificationvalues:
        if c not in classificationtypes:
            raise ValueError("%s is a value in classifications but is not in classificationtypes" % c)
    nclassifications = len(classificationtypes)
    totalxsize = float(lmargin + rmargin + xsize)
    totalysize = float(tmargin + bmargin + titlemargin + plotheight * nclassifications)
    fig = pylab.figure(figsize=(totalxsize, totalysize), facecolor='white')
    ax = pylab.axes([lmargin / totalxsize, bmargin / totalysize, 1.0 - (lmargin + rmargin) / totalxsize, 1.0 - (tmargin + bmargin + titlemargin) / totalysize])

    # plot the data
    ys = []
    xs = []
    xmin = xmax = None
    yticks = []
    for iclassification in range(nclassifications):
        ixdata = [xvalues[site] for site in sites if classifications[site] == classificationtypes[iclassification]]
        yticks.append('%s\n(%d sites)' % (classificationtypes[iclassification], len(ixdata)))
        if xmin == None:
            (xmin, xmax) = (min(ixdata), max(ixdata))
        else:
            (xmin, xmax) = (min(min(ixdata), xmin), max(max(ixdata), xmax))
        xs.append(ixdata)
        y = (iclassification + 0.5) / float(nclassifications)
        ys.append(y)
    bp = pylab.boxplot(xs, sym='b.', vert=False, positions=ys)
    # reformat box plot
    for box in bp['boxes']:
        box.set(color='blue', linewidth=1.75)
    for whisker in bp['whiskers']:
        whisker.set(color='blue', linewidth=1.75, linestyle='solid')
    for cap in bp['caps']:
        cap.set(color='blue', linewidth=1.75)
    for median in bp['medians']:
        median.set(color='red', linewidth=1.75)
    for flier in bp['fliers']:
        flier.set(marker='.', color='blue', alpha=0.5)
    # format plot
    if xmin < 0:
        xmin *= xdatamargin
    else:
        xmin /= xdatamargin
    if xmax > 0:
        xmax *= xdatamargin
    else:
        xmax /= xdatamargin
    if symmetrize_axis:
        xmax = max(abs(xmax), abs(xmin))
        xmin = -xmax
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([0.0, 1.0])
    if title:
        pylab.title(title, size=10)
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(ys))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(yticks))
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()
    


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
