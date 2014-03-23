"""Gets and aligns sequences for phylogenetic tree.

This script utilizes the `EMBOSS needle`_ program for the alignments,
and requires the ``mapmuts`` program.
"""


import sys
import os
import re
import Bio.SeqIO
import mapmuts.align
import mapmuts.sequtils


def GetUnique(seqs, mindiffs):
    """Gets unique or nearly unique sequences from a set.

    *seqs* is a list of sequences as *(header, sequence)* 2-tuples
    that are assumed to be aligned (they all must be the same length).

    *mindiffs* is an integer >= 0.

    This method returns a new list *uniqueseqs* which is like *seqs*
    that all sequences that are similar to another sequence by *mindiffs*
    or fewer differences are removed. The function goes progressively
    through *seqs* in order, and as it gets a new sequence it sees if
    it has less than or equal to *mindiffs* differences with a preceding
    sequence -- if it does, then it is discarded.
    """
    uniqueseqs = []
    for (head, seq) in seqs:
        seq = seq.upper()
        for (head2, seq2) in uniqueseqs:
            assert len(seq) == len(seq2), "Sequences differ in length. Are the aligned?"
            ndiffs = len([i for i in range(len(seq)) if seq[i] != seq2[i]])
            if ndiffs <= mindiffs:
                break
        else:
            uniqueseqs.append((head, seq))
    return uniqueseqs


def PairwiseStatistics(s1, s2):
    """Returns alignment statistics for pairwise aligned sequences.

    *s1* and *s2* are strings of the same length corresponding to two
    aligned sequences.

    The returned variable is the 2-tuple *(f_ident, f_gaps)*. *f_ident*
    is the fraction of identities among the aligned non-gap positions.
    *f_gaps* is the fraction of position in the alignment that are gaps
    in one of the sequences.
    """
    assert len(s1) == len(s2)
    n_ident = n_gaps = 0
    for (x1, x2) in zip(s1, s2):
        if x1 == '-' or x2 == '-':
            n_gaps += 1
        elif x1 == x2:
            n_ident += 1
    return (n_ident / float(len(s1)), n_gaps / float(len(s1)))


def NeedleCDSandProtAlignments(refseq, seqs, needlecmd, tempfile='_alignments.temp'):
    """Uses EMBOSS needle to align sequences in *seqs* to *refseq*.

    The sequences in *seqs* and *refseq* should be coding DNA sequences.
    They must be translatable, but ambiguous nucleotides and truncation
    of incomplete sequences are allowed.

    *refseq* is a string giving the reference sequence.

    *seqs* is a list of *(header, sequence)* 2-tuples.

    *needlecmd* is the path to the EMBOSS needle program.

    *tempfile* is the name of a file temporarily created to store
    the alignments, and then deleted.

    Returns the 2-tuple *(prot_alignment, cds_alignment)*.

    *prot_alignment* contains each of the proteins encoded in 
    *seqs* as a *(header, protsequence)* 2-tuple, with the
    protein sequence aligned to that in *refseq* and with all
    gaps relative to *resfeq* stripped away.

    *cds_alignment* is an alignment of the coding DNA sequences
    in *prot_alignment*, with the nucleotide alignments done
    according to the protein alignments.
    """
    prots = []
    heads = []
    refprot = mapmuts.sequtils.Translate([('head', refseq)])[0][1]
    for (head, seq) in seqs:
        try:
            (head, prot) = mapmuts.sequtils.Translate([(head, seq)], readthrough_n=True, truncate_incomplete=True)[0]
            heads.append(head)
            prots.append(prot)
        except:
            sys.stderr.write("PROBLEM translating sequence %s" % head)
            raise
    try:
        mapmuts.align.Needle(refprot, prots, needlecmd, 'protein', tempfile)
        alignments = mapmuts.sequtils.ReadFASTA(tempfile)
    finally:
        if os.path.isfile(tempfile):
            os.remove(tempfile)
    assert len(alignments) == 2 * len(prots) == 2 * len(heads) == 2 * len(seqs)
    prot_alignment = []
    cds_alignment = []
    for i in range(len(prots)):
        prot = prots[i]
        head = heads[i]
        seq = seqs[i][1]
        assert seqs[i][0] == head
        (refa, prota) = (alignments[2 * i][1], alignments[2 * i + 1][1])
        assert len(refa) == len(prota)
        iref = iprot = 0
        alignedprot = []
        alignedcds = []
        for (aa_ref, aa_prot) in zip(refa, prota):
            assert (aa_ref == '-' or aa_ref == refprot[iref])
            assert (aa_prot == '-' or aa_prot == prot[iprot])
            if aa_ref == '-' and aa_prot != '-':
                iprot += 1
            elif aa_prot == '-' and aa_ref != '-':
                alignedprot.append(aa_prot)
                alignedcds.append('---')
                iref += 1
            elif aa_ref != '-' and aa_prot != '-':
                alignedprot.append(aa_prot)
                alignedcds.append(seq[3 * iprot : 3 * iprot + 3])
                iref += 1
                iprot += 1
            else:
                raise ValueError("Both prots in alignment have gap")
        alignedprot = ''.join(alignedprot)
        alignedcds = ''.join(alignedcds)
        assert alignedprot == mapmuts.sequtils.Translate([(head, alignedcds)], readthrough_n=True, truncate_incomplete=True, translate_gaps=True)[0][1]
        assert len(alignedprot) == len(refprot)
        prot_alignment.append((head, alignedprot))
        cds_alignment.append((head, alignedcds))
    assert len(prot_alignment) == len(cds_alignment)
    return (cds_alignment, prot_alignment)


def main():
    """Main body of script."""

    # input / output files and script parameters
    needlecmd = 'needle' # path to EMBOSS needle, this assumes it is in search path
    refseqfile = 'TEM1_cds.fasta' # sequence to which we align
    infiles = ['Lahey_TEM_GenbankSequences.gb', 'Lahey_SHV_GenbankSequences.gb'] # input files with all sequences
    outfile = 'aligned_lactamases.fasta' # output file
    min_identity = 0.6 # minimum allowed identity
    max_gaps = 0.2 # maximum allowed gap fraction
    mindiffs = 3 # only retain sequences with more than this many differences from others


    # get reference sequence
    if not os.path.isfile(refseqfile):
        raise IOError("Could not find refseqfile of %s" % refseqfile)
    refseq = mapmuts.sequtils.ReadFASTA(refseqfile)
    if len(refseq) != 1:
        raise IOError("Failed to read exactly one sequence from refseqfile %s" % refseqfile)
    refseq = refseq[0]
    refprot = mapmuts.sequtils.Translate([refseq])[0]
    print "\nUsing a reference sequence of %s.\nThis sequence is %d nucleotides in length, and encodes a protein of %d residues." % (refseq[0], len(refseq[1]), len(refprot[1]))

    # get CDS from infiles
    cds_list = []
    for infile in infiles:
        print "\nParsing CDS sequences from %s..." % infile
        if not os.path.isfile(infile):
            raise IOError("Could not find infile %s" % infile)
        ntot = nparsed = 0
        for gb_record in Bio.SeqIO.parse(open(infile), "genbank"):
            ntot += 1
            cds = [feature for feature in gb_record.features if feature.type == 'CDS']
            if len(cds) == 1: # has exactly one CDS
                cds = cds[0]
                assert cds.strand == 1, "Currently only parses plus strand CDS"
                (startindex, endindex) = (cds.location.nofuzzy_start, cds.location.nofuzzy_end)
                assert startindex < endindex, "Starting index exceeds ending index"
                head = gb_record.description
                head = head.replace(',', '')
                head = head.replace(':', '')
                head = head.replace(' ', '_')
                head = head.replace('(', '')
                head = head.replace(')', '')
                head = 'SEQ%d_%s' % (nparsed + 1, head)
                head = head[ : 99] # make maximum header length of 100 characters
                cds = (head, str(gb_record.seq[startindex : endindex]).upper())
                try:
                    prot = mapmuts.sequtils.Translate([cds])
                    cds_list.append(cds)
                    nparsed += 1
                except:
                    pass # didn't translate
        print "Read a total of %d sequences from %s, and successfully parsed exactly one translatable CDS from %d of them." % (ntot, infile, nparsed)
    print "Read a total of %d translatable CDS sequences from all input files." %  len(cds_list)

    # purge ambiguous
    print "\nPurging any sequences with ambiguous nucleotide identities."
    m = re.compile('^[ATCGatcg]+$')
    cleanseqs = []
    for (head, seq) in cds_list:
        if m.search(seq):
            cleanseqs.append((head, seq))
    cds_list = cleanseqs
    print "Retained %d sequences after purging those with ambiguous nucleotides." % len(cds_list)

    # align to other sequences
    print "\nAligning CDS sequences via proteins..."
    (alignedseqs, prot_alignment) = NeedleCDSandProtAlignments(refseq[1], [refseq] + cds_list, needlecmd)
    print "After adding the reference sequence, there are %d aligned sequences." % len(alignedseqs)
    alignedseqs = GetUnique(alignedseqs, mindiffs=mindiffs)
    print "\nAfter retaining just aligned sequences that differ from other sequences by more than %d mutations, have %d left." % (mindiffs, len(alignedseqs))
    retained = []
    for seq in alignedseqs:
        (idents, gaps) = PairwiseStatistics(refseq[1], seq[1])
        if idents > min_identity and gaps < max_gaps:
            retained.append(seq)
    alignedseqs = retained
    print "Retained %d sequences after removing those with more than %.2f gaps or less than %.2f identities." % (len(alignedseqs), max_gaps, min_identity)
    print "Writing these aligned sequences to %s." % outfile
    mapmuts.sequtils.WriteFASTA(alignedseqs, outfile)


if __name__ == '__main__':
    main() # run the program
