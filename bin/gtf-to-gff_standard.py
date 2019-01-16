#!/usr/bin/env python

import sys

inFile = open(sys.argv[1],'r')

for line in inFile:
  #skip comment lines that start with the '#' character
  if line[0] != '#':
    #split line into columns by tab
    data = line.strip().split('\t')

    #parse the transcript/gene ID. I suck at using regex, so I usually just do a series of splits.
    transcriptID = data[-1].split('transcript_id')[-1].split(';')[0].strip()[1:-1]
    geneID = data[-1].split('gene_id')[-1].split(';')[0].strip()[1:-1]

    #replace the last column with a GFF formatted attributes columns
    #I added a GID attribute just to conserve all the GTF data
    data[-1] = "ID=" + transcriptID + ";GID=" + geneID

    #print out this new GFF line
    print '\t'.join(data)
