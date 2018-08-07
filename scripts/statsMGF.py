#!/usr/bin/env python
import sys
from collections import defaultdict

mgf_infile = open(sys.argv[1])
perc_infile = open(sys.argv[2])
crux_infile = open(sys.argv[3])

scan2title = defaultdict(set)
scan2count = defaultdict(int)
scan2demixed2data = dict()
numChimera = 0
data = set()

for line in mgf_infile:
    if "TITLE" in line:
        if len(data) > 0:
            if not scan2demixed2data.has_key(scan):
                scan2demixed2data[scan] = dict()
            scan2demixed2data[scan][demixed] = data
            data = set()

        pl = line.strip().strip('"').split(" ")
        scan = pl[0].split("=")[-1]
        demixed = pl[1].split("=")[-1]
        
        scan2count[scan]+=1
        if demixed != "original,":
            numChimera+=1

            maxIso = pl[-1].split("=")[-1]
            minIso = pl[-2].split("=")[-1]
            z = pl[-5].split("=")[-1]
            title = z+"_"+minIso+"_"+maxIso

            scan2title[scan].add(title)

    elif line[0].isdigit():
        pl = line.strip().split(" ")
        data.add(pl[0])


numScans = len(scan2count)
numDemixed = 0
numTotalDemixed = 0
avgJaccard = 0

for scan in scan2count:
    if scan2count[scan] > 1:
        numDemixed += 1
        numTotalDemixed += scan2count[scan]
        
        union = set()
        intersection = 1
        for demixed in scan2demixed2data[scan]:
            union = union.union(scan2demixed2data[scan][demixed])
            if intersection == 1:
                intersection = scan2demixed2data[scan][demixed]
            else:
                intersection = intersection.intersection(scan2demixed2data[scan][demixed])
        avgJaccard+= float(len(intersection))/len(union);

avgJaccard /= numDemixed

print numScans, numChimera, numDemixed, numTotalDemixed, avgJaccard


scan2peptide = dict()
scan2demixed2peptide = dict()
uniqueDemixed = set()
uniqueCrux = set()

perc_infile.readline()
for line in perc_infile:
    pl = line.strip().split("\t")
    scanHeader = pl[1]
    demixed = scanHeader[-2:]
    scan = scanHeader[:-2]
    peptide = pl[10]
    qval = float(pl[7])
    charge = pl[2]
    peptideC = peptide+"-"+charge

    if qval > 0.1:
        continue

    uniqueDemixed.add(peptide)

    if not scan2demixed2peptide.has_key(scan):
        scan2demixed2peptide[scan] = dict()
    scan2demixed2peptide[scan][demixed] = peptideC

crux_infile.readline()
for line in crux_infile:
    pl = line.strip().split("\t")
    scan = pl[1]
    peptide = pl[10]
    qval = float(pl[7])
    charge = pl[2]
    peptideC = peptide+"-"+charge
    
    if qval > 0.1:
        continue

    uniqueCrux.add(peptide)

    scan2peptide[scan] = peptideC


for scan in scan2peptide:
    out = [scan, scan2peptide[scan]]
    if scan2demixed2peptide.has_key(scan):
        out.append(scan2demixed2peptide[scan])
    if scan2title.has_key(scan):
        out.append(scan2title[scan])
    print out

num_unique_psms = 0

for scan in scan2demixed2peptide:
    if not scan2peptide.has_key(scan):
        print scan, scan2demixed2peptide[scan], scan2title[scan]
    uni_pep = set()
    for demixed in scan2demixed2peptide[scan]:
        uni_pep.add(scan2demixed2peptide[scan][demixed])
    num_unique_psms += len(uni_pep)
    print uni_pep

print len(uniqueDemixed), len(uniqueCrux), num_unique_psms

