.. _phyloExpCM_FreqsFromAlignment.py:

========================================
phyloExpCM_FreqsFromAlignment.py
========================================

This simple script gets the frequencies of amino acids at each site from a sequence
alignment. It then writes them to a format equivalent to that of the ``*_equilibriumpreferences.txt`` files produced by the `mapmuts`_ script ``mapmuts_inferpreferences.py``.

Typically you would want to run this script if you wanted to create a file giving the amino-acid frequencies in an alignment of sequences that could be analyzed with the `mapmuts`_ scripts ``mapmuts_siteprofileplots.py`` or ``mapmuts_preferencescorrelate.py``.

To run the script, create an input file of the format described below. Then run the script followed by the input file, as in::

    phyloExpCM_FreqsFromAlignment.py infile.txt


Format of the input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). 

The input file should contain the following keys:

* *alignmentfile* : the name of a file containing the aligned sequences. These sequences must all be aligned (of the same length). They should either be proteins or translatable DNA sequences depending on the value of *translateseqs*. This file should be in FASTA format.

* *translateseqs* : specifies whether we translate the sequences in *alignmentfile*. If *False* then the sequences in *alignmentfile* are assumed to already be protein sequences. If *True*, then the sequences in *alignmentfile* must be translatable DNA sequences. Gap codons (``---``) are translated to a ``-`` character, stop codons to a ``*`` character. In this latter case, the analysis is done on these translated sequences. An error is raised if the sequences cannot be translated.

* *includestop* : specifies whether we tally the frequencies of stop codons as a possible amino acid. If *True* then we include stop codons in the tally; if *False* then we do not.

* *pseudocounts* : integer number of pseudocounts added to the counts for each amino acid at each site. If you set this to zero, the actual fractions will be computed. But you might want to use a pseudocount (such as one) if you don't want to estimate zero frequencies.

* *outputfile* is the created file in form the of the ``*_equilibriumpreferences.txt`` file created by the `mapmuts`_ script ``mapmuts_inferpreferences.py``. Specifically, for each site all the occurrences of each amino acid is counted. Gaps are disregarded, and stop codons are either included or disregarded depending on the value of *includestop*. For each site *r* the script counts the fraction of sequences (among those being counted at this site) that have *a* as the amino acid. The results are written to the create file as in::

    #SITE   WT_AA   SITE_ENTROPY    PI_A    PI_C    PI_D    PI_E    PI_F    PI_G    PI_H    PI_I    PI_K    PI_L    PI_M    PI_N    PI_P    PI_Q    PI_R    PI_S    PI_T    PI_V    PI_W    PI_Y    PI_*
    1   na   2.85229 0.0113803   0.0402777   0.0153121   0.0320438   0.013312    0.00936795  0.026916    0.0192017   0.0224047   0.018926    0.554093    0.0445008   0.0138125   0.0212192   0.0131771   0.0152268   0.0237621   0.0102606   0.0369008   0.0367991   0.0211056
    2   na   0.968359    0.878703    0.00384431  0.00390501  0.00815666  0.0048199   0.00244426  0.00333017  0.00258064  0.00391181  0.00094265  0.0248087   0.00238593  0.00064321  0.00288323  0.00610788  0.0012339   0.001455    0.0290568   0.0107545   0.00492073  0.00311186
    3   na   2.93851 0.0295947   0.00304848  0.00227603  0.00668501  0.131168    0.000710502 0.0199653   0.00804841  0.000619715 0.366698    0.0132841   0.0223199   0.0048327   0.0170484   0.00875982  0.090712    0.0171888   0.0102176   0.177257    0.069395    0.000170873

  In this file, the ``SITE_ENTROPY`` column gives the site entropy taken to the base 2.

  In this file, all of the wildtype identities (``WT_AA``) are set to ``na`` because there is no wildtype in the sequence alignment.

  In this file, residues are numbered sequentially starting with one according to the aligned sequences in *alignmentfile*.

  Note that if *averagesites* is set to *True*, then the fractions for each site *r* will be identical, and will represent the average over the entire length of the sequence.

* *requiresubstring* is an optional argument. If it is not specified or set to *None* then nothing is done. Otherwise, it should be set to a substring -- in this case, we only compute the frequencies from sequences in *alignmentfile* that have headers that contain the substring specified by *requiresubstring*. For instance, if this substring is *HOST_Swine* then only sequences with headers that contain *HOST_swine* are counted. Note that if no sequences in *alignmentfile* contain headers with this substring, an error will be raised.

* *averagesites* is an optional Boolean switch. If it is not specified or set to *False*, then nothing is done. If you specify this option as *True*, then the fractions are **averaged over all sites**, and include the total number of occurrences of each amino acid at every site over all of the sequences. The pseudocounts specified by *pseudocounts* are then added to these total occurrences of amino acids at all sites. In this case, the reported fractions for each site in *outputfile* will be identical, since each is just reported as the average of all sites.

Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_FreqsFromAlignment.py
    alignmentfile cds_sequences.fasta
    translateseqs True
    includestop False
    pseudocounts 1
    outputfile aminoacid_frequencies.txt


Output files
----------------

Some summary output is printed to standard output. In addition, the file specified by *outputfile* is created with the format described above in the documentation for the *outputfile* specification.


.. include:: weblinks.txt
