
# import modules
import numpy as np
import MDAnalysis as md


def sumSignedInt ( OnOffArray ):
    #
    #  INPUT
    # =======
    # - OnOffArray: 1D list of signed integers
    #
    #  OUTPUT
    # ========
    # - OnOffSum: 1D list of added integer, by order and sign
    #
    
    # matrix shape 
    [ iOnOffArray , jOnOffArray ]=np.shape(OnOffArray)

    # PENDING : check array is 1D with iOnOffArray; otherwise error message
    
    
    if ( jOnOffArray > 1 ):
        # proceed for > 1 item
        # ---------------------
        
        # initialize values
        currSum=OnOffArray[0,0]
        prevBool=(currSum>0)
        OnOffSum=np.zeros((0,0), dtype=int)
        
        
        # addition loop
        for currVal in np.nditer(OnOffArray[0,1:]):
            currBool=(currVal>0)
            
            if ( currBool == prevBool ):
                currSum=currSum+currVal
            else:
                OnOffSum=np.append(OnOffSum,currSum)
                currSum=currVal
                
            prevBool=currBool

        # flush last sum
        OnOffSum=np.append(OnOffSum,currSum)

        # reshape
        lOnOffSum=len(OnOffSum)
        OnOffSum=OnOffSum.reshape(1,lOnOffSum)

    else :
        # array has only one item, just copy
        # -----------------------------------
        
        OnOffSum=OnOffArray
        
    # out
    return OnOffSum



##########




def flipCutOff ( OnOffSum, cutOff, ends=0 ):
    #
    #  INPUT
    # =======
    # - OnOffSum : 1D array of signed integers
    # - cutOff   : jump frames
    # - ends     : take into account end signs (value 1) or not (value 0)
    #             
    #
    #  OUTPUT
    # ========
    # - OnOffCutOff : 1D array, flipped sign for negative integers > cutOff
    #

    # flip to negative sign
    cutOff=-1*abs(cutOff)
    
    
    # keep start and end values for OnOffSum 
    if ends == 0 :
        OnOff1st=OnOffSum[:,0]
        OnOffLast=OnOffSum[:,-1]
        OnOffMid=OnOffSum[:,1:-1]    
    else:
        OnOffMid=OnOffSum

    
    # find values higher than cutOff
    ind=((OnOffMid>=cutOff) & (OnOffMid<=0))
    
    
    # flip sign by adding twice the cutoff values
    OnOffCutOff=OnOffMid+( OnOffMid*(ind*-2))

    
    # return start and end values to OnOffSum 
    if ends == 0 :
        OnOffCutOff=np.append(OnOff1st,OnOffCutOff)
        OnOffCutOff=np.append(OnOffCutOff,OnOffLast)
        OnOffCutOff=OnOffCutOff.reshape((1, len(OnOffCutOff) ))


    # out
    return OnOffCutOff




##########

def inOutVol ( mdUniverse, volume, selection ):
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
    # - OnOffArray : 1D array of in-out boolean values
    #
    
    # number of frames
    nFrames=len(mdUniverse.trajectory)
    
    # particle selection and volume
    strVolSel='( ' + volume + ' )' + ' and ' + '( ' + selection + ' )'
    
    # out array for bool values
    OnOffArray=np.zeros((1,nFrames), dtype=int)
    
    # iterate over trajectory
    counter=0;
    for ts in mdUniverse.trajectory:
        # intersect selection and volume
        currSel=mdUniverse.select_atoms(strVolSel)
        
        # number of atoms in selections
        currNumAtom=currSel.n_atoms
        
        # choose boolean
        if ( currNumAtom > 0 ):
            boolVolSel=1
        else:
            boolVolSel=-1
            
        # add boolean into temp array
        OnOffArray[0,counter]=boolVolSel
        counter=counter+1
        
    # out
    return OnOffArray


##########

def bool2rt ( inOutVol, jumpBack ):
    #
    #  INPUT
    # =======
    # - inOutVol  : 1D array of in-out boolean values
    # - jumpBack : blind frames 
    #
    #  OUTPUT
    # ========
    # - OnOffPos : 1D array with fileter residence times
    #

    # add values
    OnOffSum=sumSignedInt(inOutVol)

    # flip signs using jumBack
    OnOffCutOff=flipCutOff( OnOffSum, jumpBack )

    # add values after sign flip
    OnOffFrame=sumSignedInt(OnOffCutOff)

    # find positive residence times
    OnOffPos=OnOffFrame[OnOffFrame>0]

    # reorganize data
    lOnOffPos=len(OnOffPos)
    OnOffPos=OnOffPos.reshape(1,lOnOffPos)

    # out
    return OnOffPos


##########


def rtInVol ( mdUniverse, volume, selection, jumpBack ):
    #
    # residence time in volume,
    # useful for single trajectory
    #
    #  INPUT
    # =======
    # - mdUniverse : MDanalysis universe
    # - volume     : list of volumes
    # - selection  : list of selections
    # - jumpBack   : blind frames 
    #
    #
    #  OUTPUT
    # ========
    # - OnOffPos : 1D array with fileter residence times
    #

    OnOffArray=inOutVol( mdUniverse, volume, selection )
    OnOffPos=bool2rt( OnOffArray, jumpBack )

    return OnOffPos



# CHECK
def bool2frame( inOutVol, jumpBack ):
    #
    #  INPUT
    # =======
    # - inOutVol  : 1D array of in-out boolean values
    # - jumpBack : blind frames 
    #
    #  OUTPUT
    # ========
    # - outFrame : 2D array with start and end frames
    #

    # add values
    OnOffSum=sumSignedInt(inOutVol)

    # flip signs using jumBack
    OnOffCutOff=flipCutOff( OnOffSum, jumpBack )

    # add values after sign flip
    OnOffFrame=sumSignedInt(OnOffCutOff)

    # start and end indexes
    rangeFrame=np.cumsum(np.absolute(OnOffFrame))
    startFrame=np.insert(rangeFrame,0,0)
    startFrame=startFrame[0:-1]

    endFrame=rangeFrame-1
        
    # reorganize output values
    outFrame=np.zeros((0,2))
    
    iCount=0
    for oneVal in OnOffFrame[0] :        
        if oneVal > 0 :            
            oneRange=np.zeros((1,2))
            oneRange[0,0]=startFrame[iCount]
            oneRange[0,1]=endFrame[iCount]
            outFrame=np.vstack((outFrame,oneRange))            
        iCount+=1
            
    # out
    return outFrame



def frameInVol ( mdUniverse, volume, selection, jumpBack ):
    #
    # frame ranges in volume,
    # useful for single trajectory
    #
    #  INPUT
    # =======
    # - mdUniverse : MDanalysis universe
    # - volume     : list of volumes
    # - selection  : list of selections
    # - jumpBack   : blind frames 
    #
    #
    #  OUTPUT
    # ========
    # - OnOffPos : 1D array with fileter residence times
    #

    OnOffArray=inOutVol( mdUniverse, volume, selection )
    OnOffFrame=bool2frame( OnOffArray, jumpBack )

    return OnOffFrame

##########


def volFrame2bool ( inIndexFrameRange ):
    #
    # - inIndexFrameRange : [n,3] array
    #    1st col atom number
    #    2nd col start frame
    #    3rd col end frame
    #
    # - outT : tuple of 2 arrays
    #    1st array : atom or sel index
    #    2nd array  : OnOff array
    #
    
    # separate array
    inM=inIndexFrameRange[:,0]
    inM=inM.astype(int)
    frameM=inIndexFrameRange[:,1:]
    frameM=frameM.astype(int)
    
    
    # sort indexes and frame
    uniqueIndex=np.unique(inM)
    maxFrame=np.amax(frameM)

    
    # fill On-Off values
    outOnOff=np.zeros((0,maxFrame+1))

    for eachIndex in uniqueIndex :
        currOnOff=np.ones((1,maxFrame+1))
        currOnOff=currOnOff*(-1)
        
        frameRange=frameM[np.where(inM==eachIndex)]
        
        for eachRange in frameRange :
            startFrame=eachRange[0]
            endFrame=eachRange[1]
            currOnOff[0,startFrame:endFrame+1]=1
        
        outOnOff=np.vstack((outOnOff,currOnOff))

    
    # output
    uniqueIndex=uniqueIndex.reshape((1, len(uniqueIndex) ))
    outT=(uniqueIndex, outOnOff)
    return outT


##########

def volFrameJumBack ( inIndexFrameRange, jumpBack ):

    # transform volFrame data into OnOffArray
    allIndexOnOff=volFrame2bool( inIndexFrameRange )

    # separate input data
    allIndex=allIndexOnOff[0]
    allOnOff=allIndexOnOff[1]

    # output matrix
    outM=np.zeros((0,3))

    # resort data
    [iOnOff,jOnOff]=shape(allOnOff)
    
    iCount=0
    while iCount < iOnOff :

        # work with OnOff
        currIndex=allIndex[0,iCount]        
        currOnOff=allOnOff[iCount,:]
        currOnOff=currOnOff.reshape((1, len(currOnOff) ))

        # jump back
        currFrameRange=bool2frame( currOnOff, jumpBack )

        # index column
        [iCurr,jCurr]=np.shape(currFrameRange)
        colIndex=np.ones((iCurr,1))*currIndex

        # join index with frame range
        currIndexFrameRange=np.hstack((colIndex,currFrameRange ))        

        # add to output matrix
        outM=np.vstack((outM,currIndexFrameRange))

        iCount+=1

    return outM


##################
##################


##########

# PENDING
#
# 1.- visualization: PDB/DX 4 grid
# 2.- write/read output array/tupe/list
# 3.- RT over surfaces (thus atom) is more important than over volume
# 4.- X,Y,Z calculation for arbitrary volume shape should move to back burner.
#
#
#
# NOTE: 2D arrays preferred over 1D
#


##################



#################





