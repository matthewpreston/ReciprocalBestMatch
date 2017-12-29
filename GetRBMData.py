#!/usr/bin/env python
#
# Calculates the weights between query and reference sequences.
#
# For example, we have a query FASTA and an annotated reference FASTA. We wish
# to annotate the queries using the reference (where the queries are typically
# contigs). So we BLAST the query to the reference. However, this quick and
# simple procedure may not deal with paralogs as well. So one may perform the
# procedure reciprocally (BLAST the reference to the query) in order to suss out
# the identity of paralogs. This is somewhat of a manual procedure (and slightly
# heuristic). Thus I have written a document called "The Reciprocal Best Match
# Problem" (some .pdf) where I delve into the statistics of creating such an
# algorithm.
#
# One thing to note is that within that document, I have split up the RBM
# problem into 3 separate problems. I was able to (somewhat) solve the first,
# where I combine the E-values of multiple HSPs into one meaningful p-value. The
# second remains a conundrum to me - in this script, I will multiply the two
# p-values (perhaps take the log of it as they will most likely be very small).
# The last problem require graph theory, which my friend was able to solve using
# integer programming (really he used the Optimization Programming Language made
# by IBM the lazy *******). Anyways, thanks to this anonymous soul.
#
# This script will perform the algorithms derived from the 1st and 2nd problems,
# thus creating a data file required for the 3rd.
#
# This is Step 1 out of the 3 scripts in this folder.
#
# Written By: Matt Preston
# Written On: Aug  6, 2017
# Revised On: Never

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import six

#import json
from math import exp, factorial, log10, pow
import subprocess
import sys
from Bio.Blast import NCBIXML

MAX_DIG = sys.float_info.dig
MAX_NUM = 2**31 - 1 #sys.float_info.max <-- was too big

def ParseXML(xml):
    """Returs a dict of dicts containing a p-value created from HSP scores, i.e.
    
    {
        contig_1: {
                    ABCA4: 10^-7,
                    ACTN: 1.42*10^-3
                  },
        contig_2: {
                    TLR5: 4.20*10^-69
                  },
        ...
    }
    """
    result = {}
    handle = open(xml)
    for record in NCBIXML.parse(handle):
        if not record.alignments:
            continue
        lmbda, K, _ = record.ka_params
        m = record.query_length
        hits = {}
        for aln in record.alignments:
            scores = []
            for hsp in aln.hsps:
                scores.append(hsp.score)  # scores are floats
            n = aln.length
            hits[aln.hit_def] = ScoresToPValue(scores, m, n, K, lmbda)
        result[record.query] = hits
    handle.close()
    return result
    
def ScoresToPValue(scores, m, n, K, lmbda):
    """
    Given a list of scores, convert them into p-values using:
    
           p-value = 1 - e^-E * sum_{k=0}^{|O|-1} (E^k/k!)
    where
                     E = K*m*n*e^(-lmbda*s_min)
                        |O| = number of HSPs
                      m = length of query sequence
                    n = length of reference sequence
                  k,lmbda = special database constants
                         s_min = minimum score
    """
    eValue = K * m * n * exp(-lmbda * min(scores))
    A = 0
    for k in six.moves.range(len(scores)):
        A += pow(eValue, k) / factorial(k)
    #pValue = 1 - exp(-eValue) * A
    # The above, although mathematically accurate, will not yield the correct
    # answer due to double precision for floating point arithmetic being too
    # small. For example, exp(-10**-17) == 1.0, and we're going to have the
    # exponent be on the order of 10**-100. Thus the p-value it would produce
    # would always be 0.0. Instead, I will take the Taylor expansion series of
    # e^x for 1 - Ae^-E <=> 1 - A + AE - AE^2/2 + AE^3/3! - ... . This will use
    # double precision floating point numbers, which have 11 exponential bits
    # and 53 mantissa bits. Thus the summation will stop when the mantissa is no
    # longer affected by the next addition
    pValue = 1 - A + A * eValue
    old = None
    i = 2
    while old != pValue:
        old = pValue
        pValue += (-1)**(i+1) * A * pow(eValue,i) / factorial(i)
        i += 1
    return pValue

def main(argv):
    forwardXML = argv[0]
    reverseXML = argv[1]
    output = argv[2] if len(argv) > 2 else "output.csv"
    # Retrieve the p-values from query to reference
    forwardEdges = ParseXML(forwardXML)
    # Retrieve the p-values from reference to query
    reverseEdges = ParseXML(reverseXML)
    # Get pairwise weights between queries and references
    queries = list(forwardEdges.keys())
    references = list(reverseEdges.keys())
    numQueries = len(queries)
    numReferences = len(references)
    bipartite = [0] * numQueries * numReferences # Flattened 2D array
    for i, q in enumerate(queries):
        for j, r in enumerate(references):
            # Attempt to get p-values; if unsuccessful, set to 1, although it
            # should be the BLAST cutoff (default: E-value <= 5 -> p-value <= ?)
            try:
                p1 = forwardEdges[q][r]
            except KeyError as e:
                p1 = 1.0
            try:
                p2 = reverseEdges[r][q]
            except KeyError as e:
                p2 = 1.0
            # Combine bi-directional p-values into a single weight (I take the
            # -log10 of them; if the multiplication yields 0.0, then set to
            # a _really_ big number)
            try:
                value = -log10(p1 * p2)
                if value > MAX_NUM:
                    value = MAX_NUM
            except ValueError as e:
                value = MAX_NUM
            bipartite[i*numReferences+j] = value
    # Output to .dat for OPL program (using C's printf magik)
    with open("temp.dat", 'w') as o:
        o.write("qlen = %d;\n" % numQueries)
        o.write("rlen = %d;\n" % numReferences)
        o.write("weights = [\n%s\n];\n"
                % ",\n".join(["\t%.*e" % (MAX_DIG,b) for b in bipartite]))
    # Output to .json instead because you can't choose the encoding type in OPL;
    # supposedly forces the output file to be in .json format as well instead of
    # .txt which the previous was outputting (screw OPL this doesn't work)
    #data = {
    #    "parameters" : {
    #        "qlen": numQueries,
    #        "rlen": numReferences,
    #        "weights": bipartite
    #    }
    #}
    #with open(output, 'w') as o:
    #    json.dump(data, o)
    # Setup and run OPL command
    child = subprocess.Popen(["oplrun", "RBM.mod", "temp.dat"],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    child.communicate()
    # Fetch OPL data (crappily albeit but screw IBM)
    matching = []
    with open("Matching.txt") as f:
        for line in f.readlines():
            if line[:2] == "//":
                continue
            i = 0
            while i < len(line):
                if line[i].isdigit():
                    matching.append(bool(int(line[i])))
                i += 1
    # Write results to .csv
    with open(output, 'w') as o:
        o.write("Query,Reference\n")
        for i in range(len(matching)):
            if matching[i]:
                q = i // numReferences
                r = i % numReferences
                print(i, (q+1), (r+1))
                o.write("%s,%s\n" % (queries[q], references[r]))
    print(len(queries), len(references))
    # Cleanup
    subprocess.Popen(["rm", "temp.dat", "Matching.txt"]).communicate()
    sys.exit(0)

if __name__ == "__main__":
    main(sys.argv[1:])