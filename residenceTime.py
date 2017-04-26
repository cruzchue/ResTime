
# import modules
import numpy as np
import MDAnalysis as md


##########################
##########################

def selInVol ( mdUniverse, volume, selFile ):
    #
    # selections in volume
    #
    #  INPUT
    # =======
    # - mdUniverse : MDanalysis universe
    # - volume     : one volume line
    # - selFile    : selection File containing all selection lines
    #
    #
    #  OUTPUT
    # ========
    # - OnOffPos : 1D array with frame ranges (integers)
    # - outTuple : tuple that contains outLine and outList
    #              outLine : (1D array, int) containing number lines from selFile that cross volume
    #              outList : (list) short version of selFile, contains only lines that cross volume
    #
    
    # 1.- order selection data
    # -------------------------

    # read selFile into selData
    inFile=open(selFile, 'r')
    selData=inFile.read().splitlines()
    inFile.close()
    
    # selData to selTuple
    selTuple=()
    for eachLine in selData:        
        eachLine=eachLine.split()    
        intLine=np.zeros((0), dtype=int)
    
        for eachValue in eachLine:
            eachValue=int(eachValue)
            intLine=np.append(intLine, eachValue )

        selTuple=selTuple + (intLine,)

    # join selections into 1D matrix
    selMatrix=np.zeros((0), dtype=int)
    for eachLine in selTuple:
        selMatrix=np.hstack( (selMatrix , eachLine) )

    selMatrix=np.unique(selMatrix)
    selMatrix=sorted(selMatrix, key=int)

    # string line for MDanalysis
    selLine=''
    for eachVal in selMatrix:
        eachVal=str(eachVal)
        selLine=selLine+eachVal
        selLine=selLine+' '
    

    # 2.- particles in volume
    # ------------------------

    # number of frames
    nFrames=len(mdUniverse.trajectory)
        
    # particle selection and volume for MDanalysis
    strVolSel='( ' + volume + ' )' + ' and ' + '( bynum ' + selLine + ' )'

    # find particles in volume
    selInVol=np.zeros((0), dtype=int)    
    for ts in mdUniverse.trajectory:
        currSel=mdUniverse.select_atoms(strVolSel)    
        currAtoms=currSel.indices # index number, no atom number

        numAtoms=len(currAtoms)
        if numAtoms > 0 :
            currAtoms=currAtoms+1 # now, atom number 
            selInVol=np.hstack( (selInVol , currAtoms) )
        
    # sort values
    selInVol=np.unique(selInVol)
    selInVol=np.sort(selInVol)
    selCond=len(selInVol)


    # 3.- output
    # -----------

    # pick values in volume
    outList=[]
    outLine=np.zeros((0), dtype=int)
    countLine=0

    if selCond > 0 :

        for eachLine in selTuple:
            intersectVal=np.intersect1d( eachLine, selInVol)
            lenIntersect=len(intersectVal)
        
            if lenIntersect > 0 :            
                eachLine = " ".join(str(elm) for elm in eachLine)
                outList.append(eachLine)
                outLine=np.hstack( (outLine , countLine) )

            countLine+=1

    outTuple=(outLine,outList)
        
    return outTuple
    

##########

def frameInVol ( mdUniverse, volume, selection ):
    #
    # frame ranges of selection occupying volume
    #
    #  INPUT
    # =======
    # - mdUniverse : MDanalysis universe
    # - volume     : one volume
    # - selection  : one selection
    #
    #
    #  OUTPUT
    # ========
    # - OnOffPos : 1D array with frame ranges (integers)
    #

    OnOffArray=inOutVol( mdUniverse, volume, selection )
    OnOffFrame=bool2frame( OnOffArray )

    return OnOffFrame


##########


def inOutVol ( mdUniverse, volume, selection ):
    #
    # 1D On-Off array
    # +1: selection is inside volume; -1 : selection is outside volume
    #
    #  INPUT
    # =======
    # - mdUniverse : MDanalysis universe
    # - volume     : 
    # - selection  : 
    #
    #
    #  OUTPUT
    # ========
    # - OnOffArray : array with on-off values (integers)
    #
    
    # number of frames
    nFrames=len(mdUniverse.trajectory)
    
    # particle selection and volume
    strVolSel='( ' + volume + ' )' + ' and ' + '( bynum ' + selection + ' )'
    
    # out array for bool values
    OnOffArray=np.ones((nFrames), dtype=int)
    
    # iterate over trajectory
    counter=0;
    for ts in mdUniverse.trajectory:
        # intersect selection and volume
        currSel=mdUniverse.select_atoms(strVolSel)
        
        # number of atoms in selections
        currNumAtom=currSel.n_atoms

        # add on-off into out array
        if ( currNumAtom <= 0 ):
            OnOffArray[counter]=-1

        # clean
        currSel=[]
        counter=counter+1

    # out
    return OnOffArray


##########


def bool2frame( inOutVol ):
    #
    #  reshapes positive values from 1D array to frame ranges
    #
    #  INPUT
    # =======
    # - inOutVol  : 1D array of positive-negative integers
    #
    #  OUTPUT
    # ========
    # - outFrame : 2D array with start and end frames (integers)
    #

    # add values
    OnOffSum=sumSignedInt(inOutVol)

    # start and end indexes
    rangeFrame=np.cumsum(np.absolute(OnOffSum))
    
    startFrame=np.insert(rangeFrame,0,0)    
    startFrame=startFrame[0:-1]
    endFrame=rangeFrame-1
        
    # reorganize output values
    outFrame=np.zeros((0,2), dtype=int)

    iCount=0
    for oneVal in OnOffSum :        
        if oneVal > 0 :
            oneRange=np.zeros((1,2), dtype=int)
            oneRange[0,0]=startFrame[iCount]
            oneRange[0,1]=endFrame[iCount]
            outFrame=np.vstack((outFrame,oneRange))            
        iCount+=1
            
    # out
    return outFrame


def sumSignedInt ( OnOffArray ):
    #
    # add values of the same sign
    #
    #  INPUT
    # =======
    # - OnOffArray: 1D array signed integers
    #
    #  OUTPUT
    # ========
    # - OnOffSum: 1D array of added integers
    #
    
    # length of array 
    lenOnOffArray=len(OnOffArray)
    
    
    if ( lenOnOffArray > 1 ):
        # proceed for > 1 item
        # ---------------------
        
        # initialize values
        currSum=OnOffArray[0]
        prevBool=(currSum>0)
        OnOffSum=np.zeros((0), dtype=int)
        
        
        # addition loop
        for currVal in OnOffArray[1:]:
            currBool=(currVal>0)
            
            if ( currBool == prevBool ):
                currSum=currSum+currVal
            else:
                OnOffSum=np.append(OnOffSum,currSum)
                currSum=currVal
                
            prevBool=currBool

        # flush last sum
        OnOffSum=np.append(OnOffSum,currSum)

        
    else :
        # array has only one item, just copy
        # -----------------------------------
        
        OnOffSum=OnOffArray
        
    # out
    return OnOffSum


