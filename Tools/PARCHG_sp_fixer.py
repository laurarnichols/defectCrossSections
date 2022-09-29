from __future__ import division
import linecache
import math
import numpy as np
import fortranformat as ff
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b','--bands', nargs='+', help='<Required> Set band indexes to fix.', required=True, type=int)
args = parser.parse_args()

bands = args.bands
print(f'Fixing bands {bands} ...')

for iband in bands:
    inputFile = f'PARCHG.{iband:04d}.ALLK'
    outputFile = inputFile + '.fixed'
    
    f1 = open(outputFile,'w')
    for i in range(1,9):
        l = linecache.getline(inputFile,i)
        f1.write(l)
        if i == 7:
            numAtoms = sum(map(int,l.split()))
            
    for i in range(9,9+numAtoms):
        l = linecache.getline(inputFile,i)
        f1.write(l)
    
    formatted = ff.FortranRecordWriter('10(1X,E11.5)')
    l = linecache.getline(inputFile,9+numAtoms+1)
    dims = list(map(int,l.split()))
    g = int(math.ceil((dims[0]*dims[1]*dims[2])/10))
    
    for time in range(2):
        f1.write('\n')
        f1.write(l)
    
        # first scan the charge densities and fix
        # for some reason it isn't working to do the array copy with strings in the first loop
        #chgdens = np.zeros((dims[2]*dims[1]*dims[0]),dtype=str)
        print(f'read in old chgdens from {inputFile}')
        chgdens2 = []
        i2 = 0
        for i1 in range(9+numAtoms+2+time*(g+2),9+numAtoms+g+2+time*(g+2)):
            l1 = linecache.getline(inputFile,i1)
            l1 = l1.split()
            length = len(l1)
            #print l1
            #chgdens[(i2*10):(i2*10+length)] = l1
            #np.copyto(chgdens[(i2*10):(i2*10+length)],l1)
            #print chgdens[(i2*10):(i2*10+length)]
            chgdens2.extend(l1)
            i2 = i2 + 1
        chgdens = np.array(chgdens2)
        chgdens = chgdens.reshape((dims[2],dims[1],dims[0]))
        print(f'correct old chgdens from {inputFile}')
    
        for i in range(dims[2]):
            for j in range(dims[1]):
                for k in range(dims[0]):
                    #if chgdens[i][j][k] == '-':
                    #    print i,j,k
                    if chgdens[i][j][k] == '***********':
                        chgdens[i][j][k] = str((float(chgdens[i+1][j][k])+float(chgdens[i-1][j][k])+float(chgdens[i][j+1][k])+float(chgdens[i][j-1][k])+float(chgdens[i][j][k+1])+float(chgdens[i][j][k-1]))/6)
        print(f'write out fixed chgdens to {outputFile}')
        chgdens = chgdens.astype(dtype=float)
        chgdens = chgdens.flatten()
        i2 = 0
        for il in range(9+numAtoms+2,9+numAtoms+g+1):
            l1 = chgdens[(i2):(i2+10)]
            f1.write(formatted.write(l1))
            f1.write('\n')
            i2 = i2 + 10
        l1 = chgdens[(i2*10):(i2*10+length)]
        f1.write(formatted.write(l1))
        f1.write('\n')
        i2 = i2 + 10


f1.close()
#f2.close()
