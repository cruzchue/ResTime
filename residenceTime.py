
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


def flipCutOff ( OnOffSum, cutOff ):
    #
    #  INPUT
    # =======
    # - OnOffSum : 1D array of signed integers
    # - cutOff   : jump frames
    #
    #  OUTPUT
    # ========
    # - OnOffCutOff : 1D array, flipped sign for negative integers > cutOff
    #

    # flip to negative sign
    cutOff=-1*abs(cutOff)
 
    # find values higher than cutOff
    ind=((OnOffSum>=cutOff) & (OnOffSum<=0))
    
    # flip sign by adding twice the cutoff values
    OnOffCutOff=OnOffSum+( OnOffSum*(ind*-2))

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







