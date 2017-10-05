#!/usr/bin/python3

"""
1.  Read all deletions seen in  baseline
    File = sys.arg[1], X.csv, output of count_one_lane.py
    It has 3 sections:
    - Rejected (let it be for now)
    - Substitutions (currently empty because count_one_lane.py is broken)
    - Deletions -> read this section into dictionary
      Key:value break on comma
    A two level dictionary:
    - first level is individual mutations, deduced from baseline.
    - second level is counts, here enter baseline
    dict[mutation][baseline]=value
2.  In a similar way, read in deletions from N/H/MM/L counts
    - Read line and split into key:value
    - Add to dict:
     dict[mutation][filename] = value
"""

import sys
import argparse

from collections import defaultdict

from indels.ind import readDepth, readRejected, readDeletions, adjustFreq, intensity
from indels.output import printSummary

# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (
           sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect all mutations counts from different sequencing files')
    parser.add_argument('-b', '--baseline',
                        help='CSV file, first column = all mutations seen in baseline, second column = baseline counts',
                        required=True)
    parser.add_argument('-n', '--naive', help='CSV file with counts observed in the naive library', required=True)
    parser.add_argument('-hh', '--high', help='CSV file with counts observed in the positive library', required=True)
    parser.add_argument('-m', '--medium', help='CSV file with counts observed in the medium library', required=True)
    parser.add_argument('-l', '--low', help='CSV file with counts observed in the negative library', required=True)
    parser.add_argument('-d', '--depth', help='CSV file depth by position', required=True)
    # parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()

    """
    1. Read in all counts and depth
    2. Check that counts satisfy the following:
        - depth at this position is > 10K
        - mutation is observed 5 or more times in N
    3. If satisfied, add scaled entry to preferences: divide X by N & by baseline
    """
    inputs = {"baseline": args.baseline, "N": args.naive, "H": args.high, "M": args.medium, "L": args.low}
    keys = ["baseline", "N", "H", "M", "L"]
    hml = ['H', 'M', 'L']

    rejected = defaultdict(dict)
    dels = defaultdict(dict)

    depth = readDepth(args.depth)

    for name, filelocation in inputs.items():
        rejected = readRejected(name, filelocation, rejected)
        dels = readDeletions(name, filelocation, dels)

    freq = intensity(adjustFreq(depth, dels, hml), hml)
    k = ["Pos", "Residue", "Name", "SubDel", "class", "H", "M", "L"]
    printSummary(freq, k)
