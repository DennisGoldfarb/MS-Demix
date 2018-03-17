#!/usr/bin/env python
import sys

infile = open(sys.argv[1])
writeHeader = sys.argv[2]

header = infile.readline().strip().split("\t")

numCol = len(header)

numCharges = numCol - 17

newHeader = header[0:9] + header[9:9+numCharges]

for i in xrange(numCharges+1, 6):
    newHeader.append("Charge"+str(i))

newHeader += header[9+numCharges:]

if writeHeader == "T":
    print "\t".join(newHeader)

for line in infile:
    pl = line.strip().split("\t")
    newData = pl[0:9] + pl[9:9+numCharges]
    
    for i in xrange(numCharges+1, 6):
        newData.append("0")

    newData += pl[9+numCharges:]

    print "\t".join(newData) 
    
