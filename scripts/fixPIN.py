#!/usr/bin/env python
import sys
from operator import itemgetter
from math import log

infile = open(sys.argv[1])

print infile.readline().strip() # read header

scans = []

for line in infile:
    scan = line.strip().split("\t")
    scan[8] = float(scan[8]) 
    if scan[8] > 0:
        scans.append(scan)
 
scans = sorted(scans,key=itemgetter(8),reverse=True)   
scans = sorted(scans,key=itemgetter(2))

    
curr_scan = -1
curr_xcorr1 = 0

for pl in scans:
    scan = pl[2]
    XCorr = pl[8]
    
    if scan != curr_scan:
        curr_scan = scan
        curr_xcorr1 = XCorr
        
    deltaCn = 1 - ((curr_xcorr1 - XCorr) / curr_xcorr1)
    #print curr_xcorr1, XCorr, deltaCn
    #print pl
    deltaLCn = log(deltaCn)
    
    pl[6] = str(deltaLCn)
    pl[7] = str(deltaCn)
    pl[8] = str(pl[8])
    
    print "\t".join(pl)
    
