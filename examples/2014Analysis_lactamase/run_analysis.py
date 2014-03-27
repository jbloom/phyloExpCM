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

    print "Beginning execution of script at %s.\n" % time.asctime()
    if use_existing_output:
        print "Existing output will be used when possible. Note that you want to set use_existing_output to False if you want to regenerate output, such as after changing input data or analysis settings."
    else:
        print "All existing output will be deleted, and new output regenerated."

    # separator to break sections of output
    separator = '*******************************************************************' 

    # amino-acid preferences
    preferencesfile = '%s/amino_acid_preferences.txt' % os.getcwd()

    # extract the amino acid prefrences
    print "\n%s\nExtracting amino acid preferences with get_preferences.py...\n" % separator
    if os.path.isfile(preferencesfile) and use_existing_output:
        print "The preferences already exist in %s, and will not be extracted again." % preferencesfile
    else:
        os.system('python get_preferences.py')
    if not os.path.isfile(preferencesfile):
        raise ValueError("Failed to generate expected output file %s" % preferencesfile)
    print "Completed extracting the amino acid preferences with get_preferences.py into the file %s.\n%s\n\n" % (preferencesfile, separator)

    # make the preferences logo plot
    print "\n%s\nCreating logo plot of amino-acid preferences..." % separator
    logoplot = 'lactamase_site_preferences_logoplot.jpg'
    if os.path.isfile(logoplot):
        print "Logo plot %s already exists." % logoplot
    else:
        print "Creating logo plot..."
        open('mapmuts_siteprofileplots_infile.txt', 'w').write('\n'.join([
                '#Input file for mapmuts_siteprofileplots.py',
                'sitepreferences %s' % preferencesfile,
                'outfileprefix lactamase_',
                'siterange all', 
                'dsspfile 1XPB_renumbered.dssp',
                'dsspchain None',
                'add_rsa True',
                'add_ss True',
                'nperline 58',
                'includestop False',
                'sitenumbermapping sequential_to_Ambler.csv',
                ]))
        os.system('mapmuts_siteprofileplots.py mapmuts_siteprofileplots_infile.txt')
        os.system('convert -density 250 lactamase_site_preferences_logoplot.pdf %s' % logoplot)
    if not os.path.isfile(logoplot):
        raise ValueError("Failed to create logoplot %s" % logoplot)
    print "The preferences are visually displayed in %s\n%s\n\n" % (logoplot, separator)


    # build the sequence set
    alignmentfile = '%s/aligned_lactamases.fasta' % os.getcwd()
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
    raxml_tree = "%s/RAxML_bestTree.aligned_lactamases" % raxml_dir
    print "\n%s\nUsing %s to build a quick phylogenetic tree for visual analysis of potentially anomalous sequences..." % (separator, raxml)
    print "Output will be put in %s" % raxml_dir
    if use_existing_output and os.path.isdir(raxml_dir) and os.path.isfile(raxml_tree):
        print "The existing tree of %s will be used, and no new tree will be created." % raxml_tree
    else:
        if os.path.isdir(raxml_dir):
            shutil.rmtree(raxml_dir)
        os.mkdir(raxml_dir)
        os.system('%s -w %s -n aligned_lactamases -p 1 -m GTRCAT -s %s' % (raxml, raxml_dir, alignmentfile))
    if not os.path.isfile(raxml_tree):
        raise ValueError("Failed to generated expected output tree file %s" % raxml_tree)
    print "\nThe tree has been built, and is in the file %s\nYou can now visually inspect this tree for potentially anomalous sequences.\n%s\n\n" % (raxml_tree, separator)

    # build the codonPhyML trees
    script = 'phyloExpCM_runcodonPhyML.py'
    sbatch_cpus = 12 # need to set to full number on node as codonphyml is greedy
    walltime = 5 * 24 # give five days
    print "\n%s\nUsing %s to build maximum likelihood trees..." % (separator, script)
    codonphyml_trees = { # keyed by model abbreviation, value is created tree file
            'KOSI07':'%s/KOSI07_codonPhyML_tree/codonphyml_tree.newick' % os.getcwd(),
            'GY94':'%s/GY94_codonPhyML_tree/codonphyml_tree.newick' % os.getcwd(),
            }
    codonphyml_models = { # keyed by model abbreviation, value is full model specification
            'KOSI07':'KOSI07_F_omega-gamma4',
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

    # Get the subtrees
    alignmentfiles = dict([(treemodel, alignmentfile) for treemodel in codonphyml_trees.iterkeys()]) # keyed by treemodel, value is corresponding alignmentfile
    seqs = dict(mapmuts.sequtils.ReadFASTA(alignmentfile))
    subtreenames = ['SHV', 'TEM']
    print "\n%s\nGetting subtrees for %s" % (separator, ', '.join(subtreenames))
    for (maintree, maintreefile) in codonphyml_trees.items():
        print "\nExtracting subtrees for %s from %s." % (maintree, maintreefile)
        tree = Bio.Phylo.read(maintreefile, 'newick')
        print "The overall tree has %d tip nodes." % tree.count_terminals()
        for subtreename in subtreenames:
            toroot = [x for x in subtreenames if x != subtreename][0]
            nodestoroot = [node for node in tree.find_clades(name='.*%s.*' % toroot)]
            assert nodestoroot, "Found no nodes to root"
            tree.root_with_outgroup(nodestoroot)
            print "Extracting subtree with %s sequences..." % subtreename
            subtreenodes = [node for node in tree.find_clades(name='.*%s.*' % subtreename)]
            commonancestor = tree.common_ancestor(subtreenodes)
            assert commonancestor, "Cannot find a common ancestor"
            subtree = Bio.Phylo.BaseTree.Tree.from_clade(commonancestor)
            subtreeseqs = [(node.name, seqs[node.name]) for node in subtree.get_terminals()]
            print "This subtree has %d tip nodes." % subtree.count_terminals()
            subtreefile = "%s_%s.newick" % (os.path.splitext(maintreefile)[0], subtreename)
            print "Writing this tree to %s." % subtreefile
            subtreealignmentfile = "%s/aligned_%s_%s.fasta" % (os.getcwd(), maintree, subtreename)
            print "Writing the alignment for this tree to %s" % subtreealignmentfile
            mapmuts.sequtils.WriteFASTA(subtreeseqs, subtreealignmentfile)
            Bio.Phylo.write(subtree, subtreefile, 'newick')
            codonphyml_trees['%s_%s' % (maintree, subtreename)] = subtreefile
            alignmentfiles['%s_%s' % (maintree, subtreename)] = subtreealignmentfile
    print "\n%s\n\n" % separator

    # Optimize the trees for the various substitution models
    print "\n%s\nOptimizing the trees for various substitution models" % separator
    hyphy_results = {} # keyed by (treemodel, substitutionmodel) to give hyphyoutfile
    sbatch_cpus = 4 # get four CPUs as script uses lots of memory and having more CPUs typically means less other jobs on the node
    walltime =  8 * 24 # give script eight days for optimization
    processes = []
    optimizedtreedir = "%s/codonmodel_optimized_trees/" % os.getcwd()
    if not os.path.isdir(optimizedtreedir):
        os.mkdir(optimizedtreedir)
    # First the experimentally determined substitution models
    fixationmodels = ['FracTolerated', 'HalpernBruno']
    commands = {
                'hyphypath':'HYPHYMP CPU=4',
                'aapreferences':preferencesfile,
                'mutationrates':'freeparameters',
                'scalefactor':10000.0,
                'siteslist':preferencesfile,
                'keeptempfiles':'False',
                'outfileprefix':'None',
               }
    experimentalmodels = []
    for fixationmodel in fixationmodels:
        commands['fixationmodel'] = fixationmodel
        for (treemodel, tree) in codonphyml_trees.iteritems():
            commands['fastafile'] = alignmentfiles[treemodel]
            for (randomizepreferences, randstring) in [('False', ''), ('1', 'randomized')]:
                commands['randomizepreferences'] = randomizepreferences
                subdir = '%s/Tree-%s_Model-%s%s/' % (optimizedtreedir, treemodel, fixationmodel, randstring)
                hyphyoutfile = '%s/optimizedtree_results.txt' % subdir
                hyphy_results[(treemodel, "%s%s" % (fixationmodel, randstring))] = hyphyoutfile
                experimentalmodels.append("%s%s" % (fixationmodel, randstring))
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
        commands['fastafile'] = alignmentfiles[treemodel]
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
