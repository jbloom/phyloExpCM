#!python 

"""Script that runs multiple instances of ``phyloExpCM_optimizeHyphyTree.py``.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import time
import re
import multiprocessing
import subprocess
import phyloExpCM.io



def RunAnalysis(dir, hyphypath, hyphycmdfile, hyphyoutfile, hyphytreefile, hyphydistancesfile, fastafile, treefile, siteslist, model):
    """Runs ``phyloExpCM_optimizeHyphyTree.py`` script in subdirectory.

    Runs the analysis in directory *dir* after creating that directory.
    Raises an Exception if *dir* already exists.

    Will raise an Exception if the run does not complete successfully and
    create *hyphoutfile* in *dir*.

    The arguments after *dir* specify the arguments placed in the created
    ``phyloExpCM_optimizeHyphyTree_infile.txt`` within *dir*. 

    The output log file in *dir* called ``phyloExpCM_optimizeHyphyTree_log.txt``
    contains standard output from the run.

    The output file in *dir* called ``phyloExpCM_optimizeHyphyTree_errors.txt``
    contains standard error from the run.
    """
    if os.path.isdir(dir):
        raise ValueError("dir of %s already exists")
    os.mkdir(dir)
    # make path names absolute
    if os.path.isfile(hyphypath):
        # points to specific executable, make absolute
        hyphypath = os.path.abspath(hyphypath)
    fastafile = os.path.abspath(fastafile)
    treefile = os.path.abspath(treefile)
    siteslist = os.path.abspath(siteslist)
    if model[ : 12] == 'experimental' and len(model.split()) >= 2:
        (modeltype, expcmfile) = model.split(None, 2)
        expcmfile = os.path.abspath(expcmfile)
        model = "%s %s" % (modeltype, expcmfile)
    # now change to dir and run program
    os.chdir(dir)
    infilelines = [
            '# Input file for phyloExpCM_optimizeHyphyTree.py',
            'hyphypath %s' % hyphypath,
            'hyphycmdfile %s' % hyphycmdfile,
            'hyphyoutfile %s' % hyphyoutfile,
            'hyphytreefile %s' % hyphytreefile,
            'hyphydistancesfile %s' % hyphydistancesfile,
            'fastafile %s' % fastafile,
            'treefile %s' % treefile,
            'siteslist %s' % siteslist,
            'model %s' % model,
            ]
    open('phyloExpCM_optimizeHyphyTree_infile.txt', 'w').write('\n'.join(infilelines))
    subprocess.call(['phyloExpCM_optimizeHyphyTree.py', 'phyloExpCM_optimizeHyphyTree_infile.txt'], stdout=open('phyloExpCM_optimizeHyphyTree_log.txt', 'w'), stderr=open('phyloExpCM_optimizeHyphyTree_errors.txt', 'w'))
    if not os.path.isfile(hyphyoutfile):
        raise ValueError("Failed to find the expected HYPHY output file of %s in directory %s." % (hyphyoutfile, dir))


def ParseHyphyOutfile(filename):
    """Parses *hyphyoutfile* from ``phyloExpCM_optimizeHyphyTree.py``.

    *filename* gives the name of the *hyphyoutfile* created by
    ``phyloExpCM_optimizeHyphyTree.py``. An exception will be raised
    if this file does not exist.

    Returns a dictionary with the following string keys and the
    following numeric values:

        * *log_likelihood* : number giving the log likelihood 

        * *nbranchlengths* : integer number of branch lengths optimized

        * *nparameters* : integer number of parameters optimized NOT
          including the branch lengths

        * *nsharedparameters* : integer number of parameters optimized
          that are shared among all branches (such as a global
          *omega* ratio for example).
    """
    lines = open(filename).readlines()
    mapping = {'Log likelihood':float,
               'independent parameters (includes branch lengths)':int,
               'shared parameters':int,
               'number of branch lengths':int}
    d1 = {}
    for line in lines:
        key = line.split(':')[0]
        if key in mapping:
            d1[key] = mapping[key](line.split(':')[1])
            del mapping[key]
    if mapping:
        raise ValueError("Failed to find entries in %s for the following keys:\n%s" % (filename, '\n'.join(mapping.keys())))
    d = {}
    d['log_likelihood'] = d1['Log likelihood']
    d['nbranchlengths'] = d1['number of branch lengths']
    d['nsharedparameters'] = d1['shared parameters']
    d['nparameters'] = d1['independent parameters (includes branch lengths)'] - d1['number of branch lengths']
    return d



def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_multiHyphyRuns.py..."
    
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    hyphypath = phyloExpCM.io.ParseStringValue(d, 'hyphypath')
    del d['hyphypath']
    hyphycmdfile = phyloExpCM.io.ParseStringValue(d, 'hyphycmdfile')
    del d['hyphycmdfile']
    if hyphycmdfile != os.path.basename(hyphycmdfile):
        raise IOError("hyphycmdfile should be relative, not absolute, file name")
    hyphyoutfile = phyloExpCM.io.ParseStringValue(d, 'hyphyoutfile')
    del d['hyphyoutfile']
    if hyphyoutfile != os.path.basename(hyphyoutfile):
        raise IOError("hyphyoutfile should be relative, not absolute, file name")
    hyphytreefile = phyloExpCM.io.ParseStringValue(d, 'hyphytreefile')
    del d['hyphytreefile']
    if hyphytreefile != os.path.basename(hyphytreefile):
        raise IOError("hyphytreefile should be relative, not absolute, file name")
    hyphydistancesfile = phyloExpCM.io.ParseStringValue(d, 'hyphydistancesfile')
    if hyphydistancesfile != os.path.basename(hyphydistancesfile):
        raise IOError("hyphydistancesfile should be relative, not absolute, file name")
    del d['hyphydistancesfile']
    fastafile = phyloExpCM.io.ParseStringValue(d, 'fastafile')
    del d['fastafile']
    if not os.path.isfile(fastafile):
        raise IOError("Cannot find fastafile %s" % fastafile)
    treefile = phyloExpCM.io.ParseStringValue(d, 'treefile')
    del d['treefile']
    if not os.path.isfile(treefile):
        raise IOError("Cannot find treefile %s" % treefile)
    siteslist = phyloExpCM.io.ParseStringValue(d, 'siteslist')
    del d['siteslist']
    if not os.path.isfile(siteslist):
        raise IOError("Cannot file siteslist %s" % siteslist)
    summaryfile = phyloExpCM.io.ParseStringValue(d, 'summaryfile')
    del d['summaryfile']
    if os.path.isfile(summaryfile):
        print "\nRemoving existing summaryfile %s" % summaryfile
    nmulti = phyloExpCM.io.ParseIntValue(d, 'nmulti')
    del d['nmulti']
    if nmulti < 1:
        raise ValueError("nmulti must be > 1, but read value of %d" % nmulti)
    models = [(dir, model) for (dir, model) in d.iteritems()]
    if not models:
        raise ValueError("Failed to read any models")
    for (dir, model) in models:
        if os.path.isdir(dir):
            raise IOError("Directory %s already exists" % dir)
        if os.path.isfile(dir):
            raise IOError("There is already a file with the name %s" % dir)
    print "\nThis script will run phyloExpCM_optimizeHyphyTree.py for %d different models, each in a different directory." % len(models)

    # run the models nmulti at a time
    columnlabels = ['log_likelihood', 'nbranchlengths', 'nparameters', 'nsharedparameters']
    failed = []
    print "\nBeginning the analyses. Summary statistics will be written to %s..." % summaryfile
    f = open(summaryfile, 'w')
    f.write('# model, %s\n' % ', '.join(columnlabels))
    processes = {}
    start_times = {}
    for (dir, model) in models:
        print "\nRunning the analysis in subdirectory %s for model %s start at %s." % (dir, model, time.asctime())
        args = (dir, hyphypath, hyphycmdfile, hyphyoutfile, hyphytreefile, hyphydistancesfile, fastafile, treefile, siteslist, model)
        processes[dir] = multiprocessing.Process(target=RunAnalysis, args=args)
        processes[dir].start()
        start_times[dir] = time.time()
        while len(processes) >= nmulti:
            time.sleep(1)
            for key in processes.keys():
                if not processes[key].is_alive():
                    if processes[key].exitcode:
                        print "\nERROR: failed to successfully complete the analysis in %s." % key
                        failed.append(key)
                    else:
                        cumtime = (time.time() - start_times[key]) / 60.
                        print "\nSuccessfully completed the analysis in %s at %s (after %.2f hours)." % (key, time.asctime(), cumtime)
                        sys.stdout.flush()
                        d = ParseHyphyOutfile("%s/%s" % (key, hyphyoutfile))
                        f.write("%s, %g, %d, %d, %d\n" % (key, d[columnlabels[0]], d[columnlabels[1]], d[columnlabels[2]], d[columnlabels[3]]))
                        f.flush()
                    del processes[key]
                    del start_times[key]
    while processes:
        time.sleep(1)
        for key in processes.keys():
            if not processes[key].is_alive():
                if processes[key].exitcode:
                    print "\nERROR: failed to successfully complete the analysis in %s." % key
                    failed.append(key)
                else:
                    print "\nSuccessfully completed the analysis in %s." % key
                    d = ParseHyphyOutfile("%s/%s" % (key, hyphyoutfile))
                    f.write("%s, %g, %d, %d, %d\n" % (key, d[columnlabels[0]], d[columnlabels[1]], d[columnlabels[2]], d[columnlabels[3]]))
                    f.flush()
                del processes[key]
    f.close()
    if failed:
        print "\nFailed to complete analyses for the following models:\n%s\nResults for all others are summarized in %s." % (', '.join(failed), summaryfile)
    else:
        print "\nSuccessfully completed the analyses for all models. The results are summarized in %s." % summaryfile

    # script done
    print "\nScript complete."
    sys.stdout.flush()
    time.sleep(0.1)



main() # run the script
