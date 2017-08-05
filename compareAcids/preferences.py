#!/usr/bin/python3

"""
1.  Read all deletions seen in  baseline
    File = sys.arg[1], X.csv, output of indels.py
    It has 3 sections:
    - Rejected (let it be for now)
    - Substitutions (currently empty because indels.py is broken)
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
import csv
import argparse

from collections import defaultdict

# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (
           sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

parser = argparse.ArgumentParser(description='Collect all mutations counts from different sequencing files')
parser.add_argument('-b','--baseline', help='CSV file, first column = all mutations seen in baseline, second column = baseline counts',required=True)
parser.add_argument('-n','--naive', help='CSV file with counts observed in the naive library', required=True)
parser.add_argument('-hh','--high', help='CSV file with counts observed in the positive library', required=True)
parser.add_argument('-m','--medium', help='CSV file with counts observed in the medium library', required=True)
parser.add_argument('-l','--low', help='CSV file with counts observed in the negative library', required=True)
parser.add_argument('-d','--depth', help='CSV file depth by position', required=True)
#parser.add_argument('-o','--output',help='Output file name', required=True)
args = parser.parse_args()

def readDepth(depthfile):
    with open(depthfile) as f:
        d = csv.reader(f, delimiter=',')
        h = next(d)
        depth = defaultdict(dict)
        for line in d:
            depth[int(line[0])] = {h[i]:int(line[i]) for i in [1,2,3,4]}
    return depth       


""" Counts files comes in order of: Rejected, Substitutions, Deletions"""

def readRejected(name, filelocation, rejected):
    """ Read in file and add counts to corresponding index"""
    with open(filelocation, newline='\n') as f:
        rr = csv.reader(f, delimiter=',') # gives list of strings
        for line in rr:       
            if line[0] == "Rejected":
                continue
            elif line[0] == "Substitutions":
                break
            else:
                rejected[line[0]][str(name)] = int(line[1])
    return rejected


def readDeletions(name, filelocation, dels):
    """ Read in file and add counts to corresponding index"""
    with open(filelocation, newline='\n') as f:
        rr = csv.reader(f, delimiter=',') # gives list of strings
        for line in rr:
            if "Deletions" in line[0]:
                for line in rr:
                    if line[0] == "Deletions":
                        continue
                    else:
                        dels[line[0]][str(name)] = int(line[1])

    return dels

def typesubdel(mut): # mut is a strin

    m = mut.split()
    pos = int(m.pop(0))
    name = ' '.join(m)
    
    subdel = ""
    for i in range(0,len(m),3):
        if m[i+2] == '---':
            subdel += "d"
        else:
            subdel += "s"

    return pos, name, subdel


def adjustFreq(depth, dels, hml):
    # start with three interesting fractions & empty adjusted frequencies dict
    freq = defaultdict(dict)
       
    for mutation, library in dels.items():
        pos, name, subdel = typesubdel(mutation)       
        if int(library.get("N",0)) < 4:
            continue
        elif depth[pos]["N"] < 10000:
            continue
        else:
            # we have found a well characterized mutation
            freq[pos]["Name"] = name
            freq[pos]["Residue"] = int(1 + (int(mutation.split(' ')[1])-1)/3)
            freq[pos]["SubDel"] = subdel
            for entry in hml:
                freq[pos][entry] = (library.get(entry,0)*depth[pos][entry]) / (
                                    library.get("baseline")*library.get("N")*depth[pos]["N"])
    return freq


def intensity(freq, hml):
    for pos, v in freq.items():
        t = {k : v[k] for k in hml}
        m = max(t, key=lambda x: t[x])
        #max_key = max(v, key=lambda k: v[k])
        if t[m] > 1:
            freq[pos]["class"] = m
        else:
            freq[pos]["class"] = 0
    return freq     


def printSummary(dictionary, keys):
    print(','.join(keys))
    for mutation, counts in sorted(dictionary.items()):
        s = [str(counts.get(library)) for library in keys[1:]]
        print(mutation, ','.join(s), sep=",")
        
            
def main():
    """
    1. Read in all counts and depth
    2. Check that counts satisfy the following:
        - depth at this position is > 10K
        - mutation is observed 5 or more times in N
    3. If satisfied, add scaled entry to preferences: divide X by N & by baseline
    """
    inputs = {"baseline":args.baseline, "N":args.naive, "H":args.high, "M":args.medium, "L":args.low}
    keys = ["baseline", "N", "H", "M", "L"]
    hml = ['H', 'M', 'L']
    
    rejected = defaultdict(dict)
    dels = defaultdict(dict)
    
    depth = readDepth(args.depth)

    for name, filelocation in inputs.items():
        rejected = readRejected(name, filelocation, rejected)
        dels = readDeletions(name, filelocation, dels)
    
    freq = intensity(adjustFreq(depth, dels, hml),hml)
    k = ["Pos","Residue","Name","SubDel","class","H","M","L"]
    printSummary(freq,k)
      





main()
