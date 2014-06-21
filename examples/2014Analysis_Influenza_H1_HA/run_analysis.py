"""Master Python script that runs the overall analysis in this directory.

Uses ``sbatch`` and the Python ``multiprocessing`` module to run some of
the analyses in parallel for faster completion.

The locations and names of various files are hard-coded into the script.

Because this is my personal script for running jobs, the code is not
fully documented elsewhere, and you will have to look at the script source
code to understand exactly how it works. It has been tested 
on the FHCRC's computing core. 

Written by Jesse Bloom.
"""


import os
import re
import string
import time
import copy
import shutil
import glob
import multiprocessing
import Bio.Phylo
import mapmuts.sequtils
import phyloExpCM.hyphy


def RunScript(rundir, run_name, script_name, commands, use_sbatch, sbatch_cpus, walltime=None):
    """Runs a ``mapmuts`` script.

    *rundir* is the directory in which we run the job. Created if it does
    not exist.

    *run_name* is the name of the run, which should be a string without
    spaces. The input file has this prefix followed by ``_infile.txt``.

    *script_name* is the name of the script that we run.

    *commands* contains the commands written to the input file. It is a list 
    of 2-tuples giving the key / value pairs.
    Both keys and values should be strings.

    *use_sbatch* is a Boolean switch specifying whether we use ``sbatch``
    to run the script. If *False*, the script is just run with the command
    line instruction. If *True*, then ``sbatch`` is used, and the command file
    has the prefix *run_name* followed by the suffix ``.sbatch``.

    *sbatch_cpus* is an option that is only meaningful if *use_sbatch* is 
    *True*. It gives the integer number of CPUs that are claimed via
    ``sbatch`` using the option ``sbatch -c``. 

    *waltime* is an option that is only meaningful if *use_sbatch* is
    *True*. If so, it should be an integer giving the number of hours 
    to allocate for the job. If *walltime* has its default value of 
    *None*, no wall time for the job is specified.

    It is assumed that the script can be run at the command line using::

        script_name infile

    Returns *runfailed*: *True* if run failed, and *False* otherwise.
    """
    print "Running %s for %s in directory %s..." % (script_name, run_name, rundir)
    currdir = os.getcwd()
    if not os.path.isdir(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)
    if (not run_name) or not all([x not in string.whitespace for x in run_name]):
        raise ValueError("Invalid run_name of %s" % run_name)
    infile = '%s_infile.txt' % run_name
    open(infile, 'w').write('# input file for running script %s for %s\n%s' % (script_name, run_name, '\n'.join(['%s %s' % (key, value) for (key, value) in commands])))
    if use_sbatch:
        sbatchfile = '%s.sbatch' % run_name # sbatch command file
        jobidfile = 'sbatch_%s_jobid' % run_name # holds sbatch job id
        jobstatusfile = 'sbatch_%s_jobstatus' % run_name # holds sbatch job status
        joberrorsfile = 'sbatch_%s_errors' % run_name # holds sbatch job errors
        sbatch_f = open(sbatchfile, 'w')
        sbatch_f.write('#!/bin/sh\n#SBATCH\n')
        if walltime:
            sbatch_f.write('#PBS -l walltime=%d:00:00\n' % walltime)
        sbatch_f.write('%s %s' % (script_name, infile))
        sbatch_f.close()
        os.system('sbatch -c %d -e %s %s > %s' % (sbatch_cpus, joberrorsfile, sbatchfile, jobidfile))
        time.sleep(1) # short 1 second delay
        jobid = int(open(jobidfile).read().split()[-1])
        while True:
            time.sleep(1) # delay 1 second
            os.system('squeue -j %d > %s' % (jobid, jobstatusfile))
            lines = open(jobstatusfile).readlines()
            if len(lines) < 2:
                break # no longer in slurm queue
        errors = open(joberrorsfile).read().strip()
    else:
        errors = os.system('%s %s' % (script_name, infile))
    os.chdir(currdir)
    if errors:
        print "ERROR running %s for %s in directory %s." % (script_name, run_name, rundir)
        return True
    else:
        print "Successfully completed running %s for %s in directory %s." % (script_name, run_name, rundir)
        return False


def RunProcesses(processes, nmultiruns):
    """Runs a list *multiprocessing.Process* processes.

    *processes* is a list of *multiprocessing.Process* objects that
    have not yet been started. If an empty list, then just returns 
    with no action taken.

    *nmultiruns* is an integer >= 1 indicating the number of simultaneous
    processes to run.

    Runs the processes in *processes*, making sure to never have more than
    *nmultiruns* running at a time. If any of the processes fail (return
    an exitcode with a boolean value other than *False*), an exception
    is raised immediately. Otherwise, this function finishes when all
    processes have completed.
    """
    if not processes:
        return
    if not (nmultiruns >= 1 and isinstance(nmultiruns, int)):
        raise ValueError("nmultiruns must be an integer >= 1")
    processes_started = [False] * len(processes)
    processes_running = [False] * len(processes)
    processes_finished = [False] * len(processes)
    while not all(processes_finished):
        if (processes_running.count(True) < nmultiruns) and not all(processes_started):
            i = processes_started.index(False)
            processes[i].start()
            processes_started[i] = True
            processes_running[i] = True
        for i in range(len(processes)):
            if processes_running[i]:
                if not processes[i].is_alive():
                    processes_running[i] = False
                    processes_finished[i] = True
                    if processes[i].exitcode:
                        raise IOError("One of the processes failed to complete.")
        time.sleep(1)
        


def main():
    """Main body of script."""

    # If True we don't overwrite existing output; if False we regenerate everything
    use_existing_output = True

    # Do we use sbatch to submit some of the jobs?
    use_sbatch = True

    # amino acid preferences
    preferencesets = ['combined', 'replicate_1', 'replicate_2', 'replicate_3']
    preferencesfiles = ['%s/%s_amino_acid_preferences.txt' % (os.getcwd(), preferenceset) for preferenceset in preferencesets]

    # mutation rates
    mutationrates = '%s/mutationrates.txt' % os.getcwd()

    # site range
    siterange = (2, 565) # include these residues in analyses

    # List of sequence sets
    sequencesets = ['H1_HumanSwine']

    print "Beginning execution of script at %s.\n" % time.asctime()
    if use_existing_output:
        print "Existing output will be used when possible. Note that you want to set use_existing_output to False if you want to regenerate output, such as after changing input data or analysis settings."
    else:
        print "All existing output will be deleted, and new output regenerated."

    # separator to break sections of output
    separator = '*******************************************************************' 

    # build the alignments
    print "\n%s\nBuilding the sequence sets...\n" % separator
    alignmentfiles = {}
    for sequenceset in sequencesets:
        alignmentfiles[sequenceset] = '%s/%s_alignment.fasta' % (os.getcwd(), sequenceset)
        if os.path.isfile(alignmentfiles[sequenceset]) and use_existing_output:
            print "The existing alignment file %s will be used for sequence set %s." % (alignmentfiles[sequenceset], sequenceset)
        else:
            print "Creating alignment file %s for sequence set %s." % (alignmentfiles[sequenceset], sequenceset)
            os.system('python parse_%s.py' % sequenceset)
        if not os.path.isfile(alignmentfiles[sequenceset]):
            raise ValueError("Failed to generate expected output file %s" % alignmentfiles[sequenceset])
    print "\nCompleted building of sequence sets.\n%s\n\n" % separator

    # build the RAxML tree
    raxml = 'raxmlHPC-SSE3'
    raxml_dir = '%s/RAxML_output/' % os.getcwd()
    if not os.path.isdir(raxml_dir):
        os.mkdir(raxml_dir)
    print "\n%s\nUsing %s to build a quick phylogenetic tree for visual analysis of potentially anomalous sequences..." % (separator, raxml)
    for sequenceset in sequencesets:
        print "\nFor sequence set %s" % sequenceset
        raxml_tree = "%s/RAxML_bestTree.%s" % (raxml_dir, sequenceset)
        if use_existing_output and os.path.isfile(raxml_tree):
            print "The existing tree of %s will be used, and no new tree will be created." % raxml_tree
        else:
            oldfiles = glob.glob("%s/*%s*" % (raxml_dir, sequenceset))
            for fname in oldfiles:
                os.remove(fname)
            os.system('%s -w %s -n %s -p 1 -m GTRCAT -s %s' % (raxml, raxml_dir, sequenceset, alignmentfiles[sequenceset]))
        if not os.path.isfile(raxml_tree):
            raise ValueError("Failed to generated expected output tree file %s" % raxml_tree)
        print "The tree for %s has been built, and is in the file %s\nYou can now visually inspect this tree for potentially anomalous sequences.\n." % (sequenceset, alignmentfiles[sequenceset])
    print "\nCompleted RAxML analyses.\n%s\n\n" % separator

    # build the codonPhyML trees
    script = 'phyloExpCM_runcodonPhyML.py'
    sbatch_cpus = 12 # need to set to full number on node as codonphyml is greedy
    walltime = 5 * 24 # give five days
    print "\n%s\nUsing %s to build maximum likelihood trees..." % (separator, script)
    codonphyml_models = { # keyed by model abbreviation, value is full model specification
            'KOSI07':'KOSI07_F_omega-gamma4',
            'GY94':'GY94_CF3x4_omega-gamma4', 
            }
    codonphyml_trees = {}
    alignmentfiles_by_tree = {}
    processes = []
    for (model, modelspecs) in codonphyml_models.iteritems():
        for sequenceset in sequencesets:
            tree = "%s/CodonPhyML_Tree_%s_%s/codonphyml_tree.newick" % (os.getcwd(), sequenceset, model)
            codonphyml_trees['%s_%s' % (sequenceset, model)] = tree
            alignmentfiles_by_tree['%s_%s' % (sequenceset, model)] = alignmentfiles[sequenceset]
            if use_existing_output and os.path.isfile(tree):
                print "The output tree %s already exists, so no new tree will be created." % tree
            else:
                print "\nUsing codonPhyML to build the tree %s" % tree
                commands = [('seqfile', alignmentfiles[sequenceset]),
                            ('codonPhyML_path', 'codonphyml'),
                            ('seed', '1'),
                            ('outprefix', 'codonphyml'),
                            ('model', modelspecs),
                           ]
                subdirectory = os.path.dirname(tree)
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdirectory, "phyloExpCM_runcodonphyml", script, commands, use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
    RunProcesses(processes, nmultiruns=len(processes))
    for (model, tree) in codonphyml_trees.iteritems():
        if not os.path.isfile(tree):
            raise ValueError("Failed to find the expected output tree %s" % tree)
    print "\n%s\n\n" % separator

    # compute evolutionary equilibrium frequencies, make logo plots, compute correlations
    print "\n%s\nComputing evolutionary equilibrium frequencies..." % separator
    evolfreqsdir = '%s/evolutionary_frequencies/' % os.getcwd()
    if not os.path.isdir(evolfreqsdir):
        os.mkdir(evolfreqsdir)
    for (preferenceset, preferencesfile) in zip(preferencesets, preferencesfiles):
        tempfile = "%s/temp.ibf" % evolfreqsdir # create and delete this file
        evolfreqsfile = '%s/%s_evolutionary_equilibriumfreqs.txt' % (evolfreqsdir, preferenceset)
        if os.path.isfile(evolfreqsfile):
            "File %s already exists." % evolfreqsfile
        else:
            commands = [('aapreferences', preferencesfile),
                    ('mutspectrum', mutationrates),
                    ('scalefactor', '10000.0'),
                    ('makereversible', 'False'),
                    ('model', 'HalpernBruno'), # it does not matter which fixation model we choose, as they have the same stationary state
                    ('evolequilfreqs', evolfreqsfile),
                    ('hyphyExpCMs', tempfile)
                   ] 
            RunScript(evolfreqsdir, 'phyloExpCM_buildHyphyExpCM_%s' % preferenceset, 'phyloExpCM_buildHyphyExpCM.py', commands, False, 1)
            if not os.path.isfile(evolfreqsfile):
                raise ValueError("Failed to create %s" % evolfreqsfile)
            else:
                print "Successfully created %s." % evolfreqsfile
            if os.path.isfile(tempfile):
                os.remove(tempfile)
        structure_dir = '%s/PDB_structure/' % os.getcwd()
        for (numbering, numberingfile) in [('sequential_numbering', 'None'), ('H3_numbering', '%s/sequential_to_H3.txt' % structure_dir)]:
            logofile = '%s/%s_%s_evolequilfreqs_site_preferences_logoplot.jpg' % (evolfreqsdir, numbering, preferenceset)
            if os.path.isfile(logofile):
                print "Logo plot %s already exists." % logofile
                continue
            print "Making logo plot %s..." % logofile
            commands = [
                    ('outfileprefix', "%s_%s_evolequilfreqs_" % (numbering, preferenceset)),
                    ('sitepreferences', evolfreqsfile),
                    ('siterange', '%d %d' % siterange),
                    ('dsspfile', '%s/1RVX_trimer_renumbered.dssp' % structure_dir),
                    ('dsspchain', 'A'),
                    ('add_rsa', 'True'),
                    ('add_ss', 'True'),
                    ('nperline', '81'),
                    ('includestop', 'False'),
                    ('sitenumbermapping', numberingfile),
                    ]
            RunScript(evolfreqsdir, 'mapmuts_siteprofileplots_%s_%s' % (numbering, preferenceset), 'mapmuts_siteprofileplots.py', commands, False, 1)
            if not os.path.isfile("%s.pdf" % os.path.splitext(logofile)[0]):
                raise ValueError("Failed to create logofile PDF")
            os.system('convert -density 250 %s.pdf %s' % (os.path.splitext(logofile)[0], logofile))
            if not os.path.isfile(logofile):
                raise ValueError("Failed to create logofile")
            print "Successfully created logofile %s" % logofile
    print "\nCompleted computing evolutionary equilibrium frequencies.\n%s\n\n" % separator

    # correlations between preferences, conservation
    print "\n%s\nComputing and plotting correlations between preferences and frequencies.\n" % separator
    correlationdir = '%s/frequency_correlations/' % os.getcwd()
    if not os.path.isdir(correlationdir):
        os.mkdir(correlationdir)
    for sequenceset in sequencesets:
        alignmentfreqs = "%s/%s_alignment_frequencies.txt" % (correlationdir, sequenceset)
        if not os.path.isfile(alignmentfreqs):
            print "\nGetting alignment frequencies from %s and writing to %s." % (alignmentfiles[sequenceset], alignmentfreqs)
            commands = [
                ('alignmentfile', alignmentfiles[sequenceset]),
                ('translateseqs', 'True'),
                ('includestop', 'False'),
                ('pseudocounts', '0'),
                ('outputfile', alignmentfreqs),
                ]
            RunScript(correlationdir, 'phyloExpCM_FreqsFromAlignment_%s' % sequenceset, 'phyloExpCM_FreqsFromAlignment.py', commands, False, 1)
    print "Determining correlations between frequencies / preferences."
    commands = [
                ('preferencesfiles', '%s %s' % (alignmentfreqs, preferencesfiles[0])),
                ('samplenames', 'natural_frequency preference'),
                ('plotdir', correlationdir),
                ('alpha', '0.1'),
                ]
    RunScript(correlationdir, 'mapmuts_preferencescorrelate', 'mapmuts_preferencescorrelate.py', commands, False, 1)
    jpg = '%s/natural_frequency_vs_preference.jpg' % correlationdir
    pdf = "%s.pdf" % os.path.splitext(jpg)[0]
    if not os.path.isfile(jpg):
        if not os.path.isfile(pdf):
            raise ValueError("Couldn't find %s" % pdf)
        os.system('convert -density 250 %s %s' % (pdf, jpg))
        assert os.path.isfile(jpg), "Couldn't find %s" % jpg
    else:
        print "Plot %s already exists." % jpg
    print "\nCompleted computing and plotting correlations.\n%s\n\n" % separator

    # Optimize the trees for the various substitution models
    fixationmodels = ['FracTolerated', 'HalpernBruno']
    print "\n%s\nOptimizing the trees for various substitution models" % separator
    hyphy_results = {} # keyed by (treemodel, substitutionmodel) to give hyphyoutfile
    sbatch_cpus = 4 # get four CPUs as script uses lots of memory and having more CPUs typically means less other jobs on the node
    walltime =  8 * 24 # give script eight days for optimization
    processes = []
    optimizedtreedir = "%s/codonmodel_optimized_trees/" % os.getcwd()
    if not os.path.isdir(optimizedtreedir):
        os.mkdir(optimizedtreedir)
    # First the experimentally determined substitution models
    commands = {
                'hyphypath':'HYPHYMP CPU=4',
                'mutationrates':mutationrates,
                'scalefactor':10000.0,
                'keeptempfiles':'False',
                'outfileprefix':'None',
               }
    experimentalmodels = []
    for fixationmodel in fixationmodels:
        commands['fixationmodel'] = fixationmodel
        for (treemodel, tree) in codonphyml_trees.iteritems():
            commands['fastafile'] = alignmentfiles_by_tree[treemodel]
            for (randomizepreferences, randstring) in [('False', ''), ('1', 'randomized')]:
                commands['randomizepreferences'] = randomizepreferences
                for (preferenceset, preferencesfile) in zip(preferencesets, preferencesfiles):
                    commands['siteslist'] = preferencesfile
                    commands['aapreferences'] = preferencesfile
                    subdir = '%s/Tree-%s_Model-%s-%s%s/' % (optimizedtreedir, treemodel, preferenceset, fixationmodel, randstring)
                    hyphyoutfile = '%s/optimizedtree_results.txt' % subdir
                    hyphy_results[(treemodel, "%s-%s%s" % (preferenceset, fixationmodel, randstring))] = hyphyoutfile
                    experimentalmodels.append("%s-%s%s" % (preferenceset, fixationmodel, randstring))
                    commands['treefile'] = tree
                    if os.path.isfile(hyphyoutfile) and use_existing_output:
                        print "The output for the HYPHY analysis of the %s tree with the %s%s fixation model already exists in %s" % (treemodel, fixationmodel, randstring, hyphyoutfile)
                    else:
                        print "Using HYPHY to analyze the %s tree with the %s%s fixation model to create output file %s" % (treemodel, fixationmodel, randstring, hyphyoutfile)
                        processes.append(multiprocessing.Process(target=RunScript,\
                            args=(subdir, 'phyloExpCM_ExpModelOptimizeHyphyTree', 'phyloExpCM_ExpModelOptimizeHyphyTree.py', list(commands.items()), use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
    # Now the non-experimentally determined substitution models
    substitutionmodels = [
                'KOSI07_F_omega-global-one_rates-one',
                'KOSI07_F_omega-global-one_rates-gamma4',
                'KOSI07_F_omega-global-gamma4_rates-one', 
                'KOSI07_F_omega-global-gamma4_rates-gamma4',
#                'KOSI07_F_omega-branchlocal-one_rates-gamma4',
                'GY94_CF3x4_omega-global-one_rates-one',
                'GY94_CF3x4_omega-global-one_rates-gamma4',
                'GY94_CF3x4_omega-global-gamma4_rates-one',
                'GY94_CF3x4_omega-global-gamma4_rates-gamma4',
#                'GY94_CF3x4_omega-branchlocal-one_rates-gamma4',
                ]
    substitutionmodels = dict([(x, x) for x in substitutionmodels])
    print "Using %s to optimize the tree branch lengths for different substitution models and compute the resulting likelihoods..." % (script)
    commands = {
                'hyphypath':'HYPHYMP CPU=4',
                'hyphycmdfile':'hyphy_cmds.bf',
                'hyphyoutfile':'hyphy_output.txt',
                'hyphytreefile':'hyphy_tree.newick',
                'hyphydistancesfile':'None',
                'siteslist':preferencesfile,
               }
    for (treemodel, tree) in codonphyml_trees.iteritems():
        commands['fastafile'] = alignmentfiles_by_tree[treemodel]
        commands['treefile'] = tree
        for (substitutionmodel, substitutionmodelspecs) in substitutionmodels.iteritems():
            subdir = '%s/Tree-%s_Model-%s/' % (optimizedtreedir, treemodel, substitutionmodel)
            hyphyoutfile = '%s/hyphy_output.txt' % subdir
            hyphy_results[(treemodel, substitutionmodel)] = hyphyoutfile
            commands['model'] = substitutionmodelspecs
            if os.path.isfile(hyphyoutfile) and use_existing_output:
                print "The output for the HYPHY analysis of the %s tree with the %s substitution model already exists in %s" % (treemodel, substitutionmodel, hyphyoutfile)
            else:
                print "Using HYPHY to analyze the %s tree with the %s substitution model to create output file %s" % (treemodel, substitutionmodel, hyphyoutfile)
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, 'phyloExpCM_optimizeHyphyTree', 'phyloExpCM_optimizeHyphyTree.py', list(commands.items()), use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
    RunProcesses(processes, nmultiruns=len(processes))
    for ((treemodel, substitutionmodel), hyphyoutfile) in hyphy_results.iteritems():
        if not os.path.isfile(hyphyoutfile):
            raise ValueError("Failed to find expected HYPHY output for the %s tree with the %s substitution model, which should be in the file %s" % (treemodel, substitutionmodel, hyphyoutfile))
    for x in experimentalmodels:
        substitutionmodels[x] = [x]
    print "\nAll HYPHY analyses have been completed."
    print "\n%s\n\n" % separator

    # Make summary output of log likelihood and number of parameters
    empiricalparameters = [ # (unique re matching model name, number empirical parameters)
            (re.compile('GY94_CF3x4'), 9),
            (re.compile('KOSI07_F'), 60),
            (re.compile('FracTolerated'), 0),
            (re.compile('HalpernBruno'), 0),
            ]
    print "\n%s\nCreating summaries of log likelihoods and parameter counts.\n" % separator
    for treemodel in codonphyml_trees:
        fname = '%s_summary.csv' % treemodel
        flatexname = '%s_summary.tex' % treemodel
        print "Writing summary for tree %s to %s and %s..." % (treemodel, fname, flatexname)
        linelist = []
        latexlinelist = []
        for substitutionmodel in substitutionmodels:
            eparameters = None
            for (m, n) in empiricalparameters:
                if m.search(substitutionmodel):
                    if eparameters != None:
                        raise ValueError("%s matches two empiricalparameters" % substitutionmodel)
                    else:
                        eparameters = n
            if eparameters == None:
                raise ValueError("empiricalparameters not specified for %s" % substitutionmodel)
            resultfile = hyphy_results[(treemodel, substitutionmodel)]
            (loglikelihood, mlparameters) = phyloExpCM.hyphy.ExtractLikelihoodAndNParameters(resultfile)
            totparameters = mlparameters + eparameters
            aic = 2.0 * totparameters - 2.0 * loglikelihood
            textstring = "%s, %%g, %g, %d, %d, %d" % (substitutionmodel, loglikelihood, totparameters, mlparameters, eparameters)
            latexsubstitutionmodel = substitutionmodel.replace('_', ' ')
            latexstring = "%s & %%.1f & %.1f & %d (%d + %d) " % (latexsubstitutionmodel, loglikelihood, totparameters, mlparameters, eparameters)
            linelist.append((aic, textstring))
            latexlinelist.append((aic, latexstring))
        linelist.sort()
        latexlinelist.sort()
        aics = [tup[0] for tup in linelist]
        daics = [aic - min(aics) for aic in aics]
        linelist = [line % daic for ((aic, line), daic) in zip(linelist, daics)]
        latexlinelist = [line % daic for ((aic, line), daic) in zip(latexlinelist, daics)]
        f = open(fname, 'w')
        flatex = open(flatexname, 'w')
        f.write('#Summary for tree %s.\n#\n#SUBSTITUTION_MODEL, dAIC, LOG_LIKELIHOOD, FREE_PARAMETERS, MAXIMUM_LIKELIHOOD_PARAMETERS, EMPIRICAL_PARAMETERS\n%s' % (treemodel, '\n'.join([line for line in linelist])))
        flatex.write('\\begin{tabular}{c|c|c|c}\nmodel & $\Delta$AIC & log likelihood & \parbox[b]{0.9in}{\center parameters (optimized + empirical)} \\\\ \hline\n%s\n\end{tabular}' % '\\\\ \n'.join([line for line in latexlinelist]))
        f.close()
        flatex.close()
    print "\nAll summaries created.\n%s\n\n" % separator

    print "\nCompleted execution of script at %s.\n" % time.asctime()




if __name__ == '__main__':
    main() # run the script
