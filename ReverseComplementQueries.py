#!/usr/bin/env python
#
# Reverse complement queries that need to be.
#
# This is Step 2 out of the 3 scripts in this folder.
#
# Written By: Matt Preston
# Written On: Aug  9, 2017
# Revised On: Never

import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML

def main(argv):
    fastaFile = argv[0]
    xmlFile = argv[1]
    matchingFile = argv[2]
    if len(argv) > 3:
        output = argv[3]
    else:
        try:
            output, ext = fastaFile.split('.',1)
            output += "_ReverseComplement." + ext
        except:
            output = fastaFile + "_ReverseComplement.fasta"
    # Load matches into memory
    matches = []
    with open(matchingFile) as f:
        header = True
        for line in f.readlines():
            if header:
                header = False
                continue
            qID, rID = line.split(',', 1)
            matches.append(qID)
    # Load sequences into memory
    seqs = {r.id: r for r in SeqIO.parse(fastaFile, "fasta")}
    # Reverse complement those that match
    with open(xmlFile) as handle:
        for record in NCBIXML.parse(handle):
            for aln in record.alignments:
                for hsp in aln.hsps:
                    # Reverse complement query sequence if:
                    # 1) It's one of the matches
                    # 2) It is reversed during BLAST's alignment or ze reference
                    #    is reversed
                    # btw I have no clue how I should space out that if clause
                    key = record.query
                    if key in seqs.keys() \
                       and ((hsp.query_start > hsp.query_end \
                                and hsp.sbjct_start < hsp.sbjct_end) \
                            or (hsp.query_start < hsp.query_end \
                                and hsp.sbjct_start > hsp.sbjct_end)):
                        seqs[key].seq = seqs[key].seq.reverse_complement()
    # Write to file
    SeqIO.write(seqs.values(), output, "fasta")
    
if __name__ == "__main__":
    main(sys.argv[1:])
