from __future__ import division
import linecache
import math
import numpy as np
import fortranformat as ff
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-b','--bands', nargs='+', help='<Required> Set band indexes to split.', required=True, type=int)
args = parser.parse_args()

bands = args.bands
print(f'Splitting bands {bands} ...')

for iband in bands:
    inputFile = f'PARCHG.{iband:04d}.ALLK'
    if os.path.isfile(inputFile + '.fixed'):
        print('Fixed file exists.')
        inputFile += '.fixed'
    print(f'Processing {inputFile} ...')
    outputFile1 = inputFile + '.up.vasp'
    outputFile2 = inputFile + '.down.vasp'
    print(f'Writing to {outputFile1} and {outputFile2} ...')
    # outputFile1 = 'PARCHG.0431.ALLK.up'
    # outputFile2 = 'PARCHG.0431.ALLK.down'

    f1 = open(outputFile1,'w')
    f2 = open(outputFile2,'w')
    for i in range(1,9):
        l = linecache.getline(inputFile,i)
        f1.write(l)
        f2.write(l)
        if i == 7:
            numAtoms = sum(map(int,l.split()))
        
    for i in range(9,9+numAtoms):
        l = linecache.getline(inputFile,i)
        f1.write(l)
        f2.write(l)
    
    f1.write('\n')
    f2.write('\n')
    
    l = linecache.getline(inputFile,9+numAtoms+1)
    dims = [int(l.split()[0]),int(l.split()[1]),int(l.split()[2])]
    g = int(math.ceil((dims[0]*dims[1]*dims[2])/10))
    

    f1.write(l)
    f2.write(l)


    formatted = ff.FortranRecordWriter('10(1X,E11.5)')
    for i1 in range(9+numAtoms+2,9+numAtoms+g+2):
        i2 = i1 + g + 2
        l1 = linecache.getline(inputFile,i1)
        l2 = linecache.getline(inputFile,i2)
        l1 = np.array(l1.split(),dtype=float)
        l2 = np.array(l2.split(),dtype=float)
        ucd = (l1+l2)/2
        dcd = (l1-l2)/2
        f1.write(formatted.write(ucd))
        f1.write('\n')
        f2.write(formatted.write(dcd))
        f2.write('\n')
    f1.close()
    f2.close()
    print(f'{inputFile} is finished.')
