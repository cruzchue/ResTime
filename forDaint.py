

import numpy as np
import MDAnalysis as md
import time

import sys
sys.path.insert(0, '../00_Scripts/')
import residenceTime08 as rt


######### INPUT VALUES ###########

volid=int(sys.argv[1])
#volid=10

dirOut='./clean01/'

psfFile='../01_Data/eastS01.psf'
dcdFile='../01_Data/eastS01.dcd'
volFile='../01_Data/eastS01.LOGIC.volData'
#selFile='../01_Data/eastS01.SHORT.selData'
selFile='../01_Data/eastS01.selData'


############## MAIN ##############

# 1.- read files
# ----------------

# load molecule
molID=md.Universe(psfFile,dcdFile)


# read volFile
inFile=open(volFile, 'r')
volData=inFile.read().splitlines()
inFile.close()

vol=volData[volid]


# read selFile
inFile=open(selFile, 'r')
selData=inFile.read().splitlines()
inFile.close()


# 2.- residence times
# --------------------

outM=rt.frameInVolManySel( molID, vol, selData )


# 3.- output
# --------------

#  add volume column
[iOut,jOut]=np.shape(outM)

if iOut == 0 :
    outM=np.ones((1,4))*-1
else :
    colVol=np.ones((iOut,1), dtype=int)
    colVol=colVol*volid
    outM=np.hstack((colVol,outM))

    
# save text file
np.savetxt(dirOut+'/vol_'+str(volid)+'.rt',outM, fmt='%i')

