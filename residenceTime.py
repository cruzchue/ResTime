
# import modules
import numpy as np
import MDAnalysis as md



def frameInVol ( mdUniverse, volume, selection ):
    #
    # frame ranges of selection occupying volume
    #
    #  INPUT
    # =======
    # - mdUniverse : MDanalysis universe
    # - volume     : list of volumes
    # - selection  : list of selections
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
    strVolSel='( ' + volume + ' )' + ' and ' + '( ' + selection + ' )'
    
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


