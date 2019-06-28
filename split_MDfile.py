#!/usr/bin/env python

# split MD frames into individual pdbs

import sys

infile = open(sys.argv[1],'r').readlines()
print('read')
filecount = 0
output = open('frame_0000.pdb','w')
for i in infile:
    line = i.split()
    if 'END' not in line[0]:
        output.write(i)
    else:
        output.close()
        filecount+=1
        print'frame_{0:04}.pdb'.format(filecount)
        output = open('frame_{0:04}.pdb'.format(filecount),'w')