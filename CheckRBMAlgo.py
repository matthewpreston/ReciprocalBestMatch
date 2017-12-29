#!/usr/bin/env python
#
# Oh sweet jesus please work
#
# This is Step 3 out of the 3 scripts in this folder.
#
# Written By: Matt Preston
# Written On: Aug  9, 2017
# Revised On: Never

from functools import wraps
import getopt
from StringIO import StringIO
import subprocess
import sys
import threading
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline
from ThreadPool import ThreadPool

USAGE = """python CheckRBMAlgo.py [options] <query.fasta> <ref.fasta> <matching.csv>

Options:
	-d DIR, --dir=FILE
		Outputs matching query and reference alignment to DIR [None]
    -o FILE, --output=FILE
        Outputs breadth of coverage and percent identities to FILE [identities.csv]
    -t INT, --threads=INT
        Number of threads to use [1]
"""

# When a MUSCLE task fails, kill the entire program. This ensures that before
# the act of suicide, stderr is not flooded with text; i.e. only one thread will
# exclaim its impending doom before it kicks the chair away from underneath
DIE_LOCK = threading.RLock()
RETURN_CODE = 0

def CheckForDeath(function):
    """When one thread dies to an error, terminate all others"""
    @wraps(function)
    def wrapper(*args, **kwargs):
        with DIE_LOCK:
            global RETURN_CODE
            if RETURN_CODE != 0:
                sys.exit(RETURN_CODE)
        return function(*args, **kwargs)
    return wrapper
    
@CheckForDeath
def GetStatsTask(seqRecord1, seqRecord2, ID, targetDict):
    """
    Calculates the breadth of coverage and identity between two sequences,
    stores it under the key 'ID' in targetDict as a tuple (breadth, identity)
    
    Parameters:
        seqRecord1(Bio.SeqRecord.SeqRecord)
            First sequence to compare
        seqRecord2(Bio.SeqRecord.SeqRecord)
            Second sequence to compare
        ID
            Key to store under
        targetDict
            Dictionary to store identities
    """
    records = [seqRecord1, seqRecord2]
    # Align with MUSCLE
    cline = MuscleCommandline()
    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    SeqIO.write(records, child.stdin, "fasta")
    stdout, stderr = child.communicate()
    # Check for errors
    with DIE_LOCK:
        # Check if other threads had errors
        global RETURN_CODE
        if RETURN_CODE:
            sys.exit(RETURN_CODE)
        # Check for local errors
        if child.returncode:
            sys.stderr.write(stderr)
            RETURN_CODE = child.returncode
            sys.exit(RETURN_CODE)
    align = AlignIO.read(StringIO(stdout), "fasta")
    # Find breadth of coverage
    breadth = float(sum(s1 != '-' and s2 != '-' for s1,s2 in zip(*align))) \
            / sum(s2 != '-' for s2 in align[1]) * 100
    # Find sequence identity between reference and consensus
    identity = float(sum(s1 == s2 for s1,s2 in zip(*align))) \
             / align.get_alignment_length() * 100
    targetDict[ID] = (breadth, identity)
    
def OutputAlignmentTask(output, *seqRecords):
    SeqIO.write(seqRecords, output, "fasta")

def main(argv):
    optstr = "hd:o:t:"
    optlist = ["help","output=","threads="]
    try:
        opts, args = getopt.getopt(argv, optstr, optlist)
    except getopt.GetoptError:
        sys.stderr.write(USAGE)
        sys.exit(-1)
    directory = None
    output = "output.csv"
    threads = 1
    for opt, arg in opts:
        if opt in ["-h", "--help"]:
            sys.stdout.write(USAGE)
            sys.exit()
        elif opt in ["-d", "--dir"]:
            directory = arg
        elif opt in ["-o", "--output"]:
            output = arg
        elif opt in ["-t", "--threads"]:
            try:
                threads = int(arg)
            except ValueError as e:
                sys.stderr.write("Invalid number of threads: %s\n" % arg)
                sys.exit(-1)
    if len(args) < 3: 
        sys.stderr.write(USAGE)
        sys.exit(-2)
        
    queryFasta = args[0]
    refFasta = args[1]
    matchingFile = args[2]
    # Make alignments directory if need be
    if directory is not None:
        subprocess.call(["mkdir", directory])
    # Set up a pool of threads for concurrent runs of MUSCLE
    pool = ThreadPool(threads-1)
    # Get queries and references from FASTAs
    queries = {r.id: r for r in SeqIO.parse(queryFasta, "fasta")}
    refs = {r.id.split('_',1)[0]: r for r in SeqIO.parse(refFasta, "fasta")}
    # Match query to reference using matching data and get breadth of coverage
    # and identity statistics
    data = {}
    with open(matchingFile) as f:
        header = True
        for line in f.readlines():
            if header:
                header = False
                continue
            qID, rID = line.split(',', 1)
            rID = rID.split('_',1)[0].strip()
            pool.AddTask(GetStatsTask,
                         queries[qID], refs[rID], rID, data)
            if directory is not None:
                pool.AddTask(OutputAlignmentTask,
                             directory + "/" + rID + ".fasta",
                             queries[qID], refs[rID])
    pool.HaveMainThreadWork()
    # Write data to .csv
    with open(output, 'w') as o:
        o.write("Gene,Breadth,Identity\n")
        for k in sorted(data.keys(), key=lambda k: k.lower()):
            b,i = data[k]
            o.write("%s,%s,%s\n" % (k,b,i))
    sys.exit(0)
    
if __name__ == "__main__":
    main(sys.argv[1:])