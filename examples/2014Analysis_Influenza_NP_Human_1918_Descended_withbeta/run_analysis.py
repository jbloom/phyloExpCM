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
import multiprocessing
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
        

def FormatModelName(name):
    """Returns a LaTex formatted model name."""
    formatted_names = {
            'fitbetaHalpernBruno':'experimental, \\ref{eq:Frxy_HalpernBruno}, free $\\beta$',
            'HalpernBruno':'experimental, \\ref{eq:Frxy_HalpernBruno}, $\\beta = 1$',
            'fitbetaFracTolerated':'experimental, \\ref{eq:Frxy_FracTolerated}, free $\\beta$',
            'FracTolerated':'experimental, \\ref{eq:Frxy_FracTolerated}, $\\beta = 1$',
            'GY94_CF3x4_omega-global-gamma4_rates-gamma4':'GY94, gamma $\omega$, gamma rates',
            'KOSI07_F_omega-global-gamma4_rates-gamma4':'KOSI07+F, gamma $\omega$, gamma rates',
            'GY94_CF3x4_omega-global-gamma4_rates-one':'GY94, gamma $\omega$, one rate',
            'KOSI07_F_omega-global-gamma4_rates-one':'KOSI07+F, gamma $\omega$, one rate',
            'GY94_CF3x4_omega-global-one_rates-gamma4':'GY94, one $\omega$, gamma rates',
            'KOSI07_F_omega-global-one_rates-gamma4':'KOSI07+F, one $\omega$, gamma rates',
            'KOSI07_F_omega-global-one_rates-one':'KOSI07+F, one $\omega$, one rate',
            'GY94_CF3x4_omega-global-one_rates-one':'GY94, one $\omega$, one rate',
            'fitbetaFracToleratedrandomized':'randomized, \\ref{eq:Frxy_FracTolerated}, free $\\beta$',
            'fitbetaHalpernBrunorandomized':'randomized, \\ref{eq:Frxy_HalpernBruno}, free $\\beta$',
            'avgaafreqsfitbeta_FracTolerated':'avg. frequencies, \\ref{eq:Frxy_FracTolerated}, free $\\beta$',
            'avgaafreqsfitbeta_HalpernBruno':'avg. frequencies, \\ref{eq:Frxy_HalpernBruno}, free $\\beta$',
            'avgaafreqs_FracTolerated':'avg. frequencies, \\ref{eq:Frxy_FracTolerated}, $\\beta = 1$',
            'avgaafreqs_HalpernBruno':'avg. frequencies, \\ref{eq:Frxy_HalpernBruno}, $\\beta = 1$',
            'FracToleratedrandomized':'randomized, \\ref{eq:Frxy_FracTolerated}',
            'HalpernBrunorandomized':'randomized, \\ref{eq:Frxy_HalpernBruno}',
                      }
    if name in formatted_names:
        return formatted_names[name]
    else:
        return name.replace('_', ' ')


def main():
    """Main body of script."""

    # If True we don't overwrite existing output; if False we regenerate everything
    use_existing_output = True

    # Do we use sbatch to submit some of the jobs?
    use_sbatch = True

    # density for JPGs created from PDFs
    jpg_density = 150

    # siterange gives the starting and ending residue analyze, inclusive
    siterange = (2, 498) # analyze residues 2 to 498, in sequential 1, 2, ... numbering

    # mutation rates
    mutationrates = '%s/mutspectrum.txt' % os.getcwd()

    # amino-acid preferences
    preferencesfile = '%s/p1_equilibriumpreferences.txt' % os.getcwd()


    print "Beginning execution of script at %s.\n" % time.asctime()
    if use_existing_output:
        print "Existing output will be used when possible. Note that you want to set use_existing_output to False if you want to regenerate output, such as after changing input data or analysis settings."
    else:
        print "All existing output will be deleted, and new output regenerated."

    # separator to break sections of output
    separator = '*******************************************************************' 

    # build the sequence set
    alignmentfile = '%s/Aligned_NPs.fasta' % os.getcwd()
    print "\n%s\nBuilding the sequence set with get_treeseqs.py...\n" % separator
    if os.path.isfile(alignmentfile) and use_existing_output:
        print "The existing alignment file of %s will be used, and no new sequence set will be created." % alignmentfile
    else:
        os.system('python get_treeseqs.py')
    if not os.path.isfile(alignmentfile):
        raise ValueError("Failed to generate expected output file %s" % alignmentfile)
    print "\nCompleted building of sequence set with get_treeseqs.py to generate the alignment file %s.\n%s\n\n" % (alignmentfile, separator)

    # build the RAxML tree
    raxml = 'raxmlHPC-SSE3'
    raxml_dir = '%s/RAxML_output/' % os.getcwd()
    raxml_tree = "%s/RAxML_bestTree.Aligned_NPs" % raxml_dir
    print "\n%s\nUsing %s to build a quick phylogenetic tree for analysis of anomalous sequences that are clear outliers from the molecular clock..." % (separator, raxml)
    print "Output will be put in %s" % raxml_dir
    if use_existing_output and os.path.isdir(raxml_dir) and os.path.isfile(raxml_tree):
        print "The existing tree of %s will be used, and no new tree will be created." % raxml_tree
    else:
        if os.path.isdir(raxml_dir):
            shutil.rmtree(raxml_dir)
        os.mkdir(raxml_dir)
        os.system('%s -w %s -n Aligned_NPs -p 1 -m GTRCAT -s %s' % (raxml, raxml_dir, alignmentfile))
    if not os.path.isfile(raxml_tree):
        raise ValueError("Failed to generated expected output tree file %s" % raxml_tree)
    print "\nThe tree has been built, and is in the file %s\nYou can now visually inspect this tree with Path-O-Gen to remove any outliers from the molecular clock.\n%s\n\n" % (raxml_tree, separator)

    # build the codonPhyML trees
    script = 'phyloExpCM_runcodonPhyML.py'
    sbatch_cpus = 12 # need to set to full number on node as codonphyml is greedy
    walltime = 4 * 24 # give four days
    print "\n%s\nUsing %s to build maximum likelihood trees..." % (separator, script)
    codonphyml_trees = { # keyed by model abbreviation, value is created tree file
            'GY94':'%s/GY94_codonPhyML_tree/codonphyml_tree.newick' % os.getcwd(),
            }
    codonphyml_models = { # keyed by model abbreviation, value is full model specification
            'GY94':'GY94_CF3x4_omega-gamma4', 
            }
    processes = []
    for (model, modelspecs) in codonphyml_models.iteritems():
        tree = codonphyml_trees[model]
        if use_existing_output and os.path.isfile(tree):
            print "The output tree %s for model %s already exists, so no new tree will be created." % (tree, model)
        else:
            print "\nUsing codonPhyML to use the model %s to build the tree %s" % (model, tree)
            commands = [('seqfile', alignmentfile),
                        ('codonPhyML_path', 'codonphyml'),
                        ('seed', '1'),
                        ('outprefix', 'codonphyml'),
                        ('model', modelspecs),
                       ]
            subdirectory = os.path.dirname(tree)
            processes.append(multiprocessing.Process(target=RunScript,\
                args=(subdirectory, "phyloExpCM_runcodonphyml", script, copy.deepcopy(commands), use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
    RunProcesses(processes, nmultiruns=len(processes))
    for (model, tree) in codonphyml_trees.iteritems():
        if not os.path.isfile(tree):
            raise ValueError("Failed to find the expected output tree %s" % tree)
    print "\n%s\n\n" % separator

    # Get the average amino-acid frequencies for each alignment
    print "\nExtracting average amino-acid frequencies for the alignment (averaged over all sites)..."
    avg_aa_freqs = '%s/avg_aa_freqs.txt' % os.getcwd()
    commands = [('alignmentfile', alignmentfile),
                ('translateseqs', 'True'),
                ('includestop', 'False'),
                ('pseudocounts', '1'),
                ('averagesites', 'True'),
                ('outputfile', avg_aa_freqs),
                ]
    if os.path.isfile(avg_aa_freqs) and use_existing_output:
        print "Using existing average amino-acid frequencies in the file %s" % avg_aa_freqs
    else:
        print "Extracting average amino-acid frequencies to %s" % avg_aa_freqs
        RunScript('./', 'phyloExpCM_FreqsFromAlignment', 'phyloExpCM_FreqsFromAlignment.py', commands, False, 1)
    if not os.path.isfile(avg_aa_freqs):
        raise ValueError("Failed to generate expected output file of %s" % avg_aa_freqs)

    # Optimize the trees for the different substitution models
    #fixationmodels = ['FracTolerated', 'HalpernBruno']
    fixationmodels = ['FracTolerated']
    print "\n%s\nOptimizing the trees for various substitution models" % separator
    hyphy_results = {} # keyed by (treemodel, substitionmodel) to give hyphyoutfile
    sbatch_cpus = 4 # get four CPUs as script uses lots of memory and having more CPUs typically means less other jobs on the node
    walltime =  5 * 24 # give four days for optimization
    processes = []
    optimizedtreedir = "%s/codonmodel_optimized_trees/" % os.getcwd()
    if not os.path.isdir(optimizedtreedir):
        os.mkdir(optimizedtreedir)
    # First the experimentally determined substitution models
    commands = {
                'hyphypath':'HYPHYMP CPU=4',
                'mutationrates':mutationrates,
                'scalefactor':10000.0,
                'siteslist':preferencesfile,
                'keeptempfiles':'False',
                'outfileprefix':'None',
                'fastafile':alignmentfile,
                'persitelikelihoods':'True',
                }
    experimentalmodels = []
    for fixationmodel in fixationmodels:
        commands['fixationmodel'] = fixationmodel
        for (treemodel, tree) in codonphyml_trees.iteritems():
            commands['treefile'] = tree
            for (fitbeta, fitbetastring) in [('False', ''), ('freeparameter', 'fitbeta')]:
                commands['fitbeta'] = fitbeta
                commands['aapreferences'] = preferencesfile
#                for (randomizepreferences, randstring) in [('False', ''), ('1', 'randomized')]:
                for (randomizepreferences, randstring) in [('False', '')]:
                    commands['randomizepreferences'] = randomizepreferences
                    istring = '%s%s%s' % (fitbetastring, fixationmodel, randstring)
                    subdir = '%s/Tree-%s_Model-%s/' % (optimizedtreedir, treemodel, istring)
                    hyphyoutfile = '%s/optimizedtree_results.txt' % subdir
                    if os.path.isfile(hyphyoutfile) and use_existing_output:
                        print "The output for the HYPHY analysis of the %s tree with the %s fixation model already exists in %s" % (treemodel, istring, hyphyoutfile)
                    else:
                        print "Using HYPHY to analyze the %s tree with the %s fixation model to create output file %s" % (treemodel, istring, hyphyoutfile)
                        processes.append(multiprocessing.Process(target=RunScript,\
                            args=(subdir, 'phyloExpCM_ExpModelOptimizeHyphyTree', 'phyloExpCM_ExpModelOptimizeHyphyTree.py', list(commands.items()), use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
                    experimentalmodels.append(istring)
                    hyphy_results[(treemodel, istring)] = hyphyoutfile
                # now add models using average aa frequencies
                commands['randomizepreferences'] = 'False'
                commands['aapreferences'] = avg_aa_freqs
                istring = 'avgaafreqs%s_%s' % (fitbetastring, fixationmodel)
                subdir = '%s/Tree-%s_Model-%s/' % (optimizedtreedir, treemodel, istring)
                hyphyoutfile = '%s/optimizedtree_results.txt' % subdir
                if os.path.isfile(hyphyoutfile) and use_existing_output:
                    print "The output for the HYPHY analysis of the %s tree with the %s fixation model already exists in %s" % (treemodel, istring, hyphyoutfile)
                else:
                    print "Using HYPHY to analyze the %s tree with the %s fixation model to create output file %s" % (tree, istring, hyphyoutfile)
                    processes.append(multiprocessing.Process(target=RunScript,\
                        args=(subdir, 'phyloExpCM_ExpModelOptimizeHyphyTree', 'phyloExpCM_ExpModelOptimizeHyphyTree.py', list(commands.items()), use_sbatch, sbatch_cpus), kwargs={'walltime':walltime}))
                hyphy_results[(treemodel, istring)] = hyphyoutfile
                experimentalmodels.append(istring)
    # Now the non-experimentally determined models
    substitutionmodels = [
                'KOSI07_F_omega-global-one_rates-one',
                'KOSI07_F_omega-global-one_rates-gamma4',
                'KOSI07_F_omega-global-gamma4_rates-one', 
                'KOSI07_F_omega-global-gamma4_rates-gamma4', 
#                'KOSI07_F_M2a_rates-one',
#                'KOSI07_F_M8_rates-one',
                'GY94_CF3x4_omega-global-one_rates-one',
                'GY94_CF3x4_omega-global-one_rates-gamma4',
                'GY94_CF3x4_omega-global-gamma4_rates-one',
                'GY94_CF3x4_omega-global-gamma4_rates-gamma4',
#                'GY94_CF3x4_M2a_rates-one',
#                'GY94_CF3x4_M8_rates-one',
                ]
    substitutionmodels = dict([(x, x) for x in substitutionmodels])
    print "\n%s\nOptimizing the tree branch lengths for different substitution models and compute the resulting likelihoods..." % separator
    commands = {
                'hyphypath':'HYPHYMP CPU=1',
                'hyphycmdfile':'hyphy_cmds.bf',
                'hyphyoutfile':'hyphy_output.txt',
                'hyphytreefile':'hyphy_tree.newick',
                'hyphydistancesfile':'None',
                'fastafile':alignmentfile,
                'siteslist':preferencesfile,
                'persitelikelihoods':'sitelikelihoods.txt',
               }
    for (treemodel, tree) in codonphyml_trees.iteritems():
        commands['treefile'] = tree
        for (substitutionmodel, substitutionmodelspecs) in substitutionmodels.iteritems():
            subdir = '%s/Tree-%s_Model-%s/' % (optimizedtreedir, treemodel, substitutionmodel)
            hyphyoutfile = '%s/hyphy_output.txt' % subdir
            if os.path.isfile(hyphyoutfile) and use_existing_output:
                print "The output for the HYPHY analysis of the %s tree with the %s substitution model already exists in %s" % (treemodel, substitutionmodel, hyphyoutfile)
            else:
                print "Using HYPHY to analyze the %s tree with the %s substitution model to create output file %s" % (treemodel, substitutionmodel, hyphyoutfile)
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, 'phyloExpCM_optimizeHyphyTree', 'phyloExpCM_optimizeHyphyTree.py', list(commands.items()), use_sbatch, 1), kwargs={'walltime':walltime}))
            hyphy_results[(treemodel, substitutionmodel)] = hyphyoutfile
            commands['model'] = substitutionmodelspecs
    RunProcesses(processes, nmultiruns=len(processes))
    for ((treemodel, substitutionmodel), hyphyoutfile) in hyphy_results.iteritems():
        if not os.path.isfile(hyphyoutfile):
            raise ValueError("Failed to find expected HYPHY output for the %s tree with the %s substitution model, which should be in the file %s" % (treemodel, substitutionmodel, hyphyoutfile))
    for x in experimentalmodels:
        substitutionmodels[x] = x
    print "\nAll HYPHY analyses have been completed."
    print "\n%s\n\n" % separator

    # Make summary output of log likelihood and number of parameters
    empiricalparameters = [ # (unique re matching model name, number empirical parameters)
            (re.compile('GY94_CF3x4'), 9),
            (re.compile('KOSI07_F'), 60),
            (re.compile('^FracTolerated'), 0),
            (re.compile('^HalpernBruno'), 0),
            (re.compile('^avgaafreqs_'), 60),
            (re.compile('^fitbeta'), 0),
            (re.compile('^avgaafreqsfitbeta_'), 60),
            ]
    print "\n%s\nCreating summaries of log likelihoods and parameter counts.\n" % separator
    parameters = { # keyed by parameter strings in HYPHY output, value is printed parameter name
                'beta':'$\\beta$' ,
                'R_AG':'$R_{A \\shortrightarrow G}$',
                'R_AT':'$R_{A \\shortrightarrow T}$',
                'R_CA':'$R_{C \\shortrightarrow A}$',
                'R_CG':'$R_{C \\shortrightarrow G}$',
                'global rate_alpha':'rate shape',
                'global omega':'$\omega$',
                'global omega_mean':'mean $\omega$',
                'global kappa':'$\kappa$',
                'global omega_alpha':'$\omega$ shape',
                }
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
            treeparameters = dict([(parameter, None) for parameter in parameters.iterkeys()])
            if 'fitbeta' not in substitutionmodel:
                del treeparameters['beta']
            phyloExpCM.hyphy.ExtractValues(resultfile, treeparameters, allowmissing=True)
            (loglikelihood, mlparameters) = phyloExpCM.hyphy.ExtractLikelihoodAndNParameters(resultfile)
            totparameters = mlparameters + eparameters
            aic = 2.0 * totparameters - 2.0 * loglikelihood
            textstring = "%s, %%g, %g, %d, %d, %d" % (substitutionmodel, loglikelihood, totparameters, mlparameters, eparameters)
            latexsubstitutionmodel = FormatModelName(substitutionmodel)
            treeparameters = ['%s = %.1f' % (parameters[parameter], value) for (parameter, value) in treeparameters.iteritems() if value]
            treeparameters.sort()
            treeparameters = ', '.join(treeparameters)
            if not treeparameters:
                treeparameters = 'none'
            latexstring = "%s & %%.1f & %.1f & %d (%d + %d) & %s" % (latexsubstitutionmodel, loglikelihood, totparameters, mlparameters, eparameters, treeparameters)
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
        flatex.write('{\\scriptsize\n\\begin{tabular}{ccccc}\nmodel & $\Delta$AIC & \parbox[b]{0.53in}{\center log likelihood} & \parbox[b]{0.7in}{\center parameters (optimized + empirical)} & optimized parameters \\\\ \hline\n%s\n\end{tabular}}' % '\\\\ \n'.join([line for line in latexlinelist]))
        f.close()
        flatex.close()
    print "\nAll summaries created.\n%s\n\n" % separator

    # compare site likelihoods
    print "\n%s\nComparing site likelihoods for the best experimental and traditional model.\n\n" % separator
    (model1name, model1file) = ('experimental', 'codonmodel_optimized_trees/Tree-GY94_Model-fitbetaFracTolerated/sitelikelihoods.txt')
    (model2name, model2file) = ('GY94', 'codonmodel_optimized_trees/Tree-GY94_Model-GY94_CF3x4_omega-global-gamma4_rates-gamma4/sitelikelihoods.txt')
    print "The two models being compared are %s (%s) versus %s (%s)" % (model1name, model1file, model2name, model2file)
    commands = [('sitelikelihoodfiles', '%s %s' % (model1file, model2file)),
                ('modelnames', '%s %s' % (model1name, model2name)),
                ('dsspfile', 'DSSP_analysis/2IQH_monomerC.dssp'),
                ('dsspchain', 'None'),
                ('outfileprefix', 'None'),
                ]
    RunScript('./', 'phyloExpCM_SiteLikelihoodComparison', 'phyloExpCM_SiteLikelihoodComparison.py', commands, False, 1)
    for x in ['SS', 'RSA']:
        os.system('convert -density 300 sitelikelihoodcomparison_by%s.pdf sitelikelihoodcomparison_by%s.jpg' % (x, x))
    print "\nCompleted the site likelihood comparison.\n%s\n\n" % separator

    print "\nCompleted execution of script at %s.\n" % time.asctime()



if __name__ == '__main__':
    main() # run the script
