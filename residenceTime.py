
# import modules
import numpy as np
import MDAnalysis as md
import time


##########################
##########################

def frameInVolManySel ( mdUniverse, volume, selectionList ):

    # 1.- pick selections in selectionList
    # -------------------------------------
    
    #'''''
    np.disp('start sorting')
    tstart=time.time()
    #'''''
    
    [ lineInSelection , inSelection ]=selInVol( mdUniverse, volume, selectionList )
    
    #'''''
    tend=time.time()
    dt=tend-tstart
    np.disp(dt)
    np.disp('finish sorting')
    #'''''


    # 2.-  matrix with index and frame ranges
    # ----------------------------------------

    outM=np.zeros((0,3), dtype=int)
    
    for lineSel,sel in zip(lineInSelection, inSelection):
        
        #'''''
        tstart=time.time()
        #'''''
    
        currFrame=frameInVol( mdUniverse , volume, sel )
        [iFrame,jFrame]=np.shape(currFrame)
        
        if iFrame > 0 :
            
            colIndex=np.ones((iFrame,1), dtype=int)
            colIndex=colIndex*lineSel
            
            currOut=np.hstack((colIndex,currFrame))
            outM=np.vstack((outM,currOut))

        #'''''
        tend=time.time()
        dt=tend-tstart
        np.disp(dt)
        #'''''


    # 3.- add volume column
    # ----------------------
    
    [iOut,jOut]=np.shape(outM)

    if iOut == 0 :
        outM=np.ones((1,3))*-1

    return outM



##########################
##########################



def selInVol ( mdUniverse, volume, selectionList ):
    #
    # selections in volume
    #
    #  INPUT
    # =======
    # - mdUniverse    : MDanalysis universe
    # - volume        : one volume line
    # - selectionList : list containing all selection lines
    #
    #
    #  OUTPUT
    # ========
    # - outTuple : tuple that contains outSelLine and outSel
    #              outSelLine : (1D array, int) containing number lines from selFile that cross volume
    #              outSel     : (list) short version of selFile, contains only lines that cross volume
    #

    #>>>>>>>
    
    # 1.- order selection data
    # -------------------------

    #>>> # read selFile into selectionList
    #>>> inFile=open(selFile, 'r')
    #>>> selectionList=inFile.read().splitlines()
    #>>> inFile.close()

    #>>>>>>>
    
    outSel=[]
    outSelLine=np.zeros((0), dtype=int)
    countLine=0


    # 2.- particles in volume
    # ------------------------

    # run over selections
    for sel in selectionList :
    
        strSel='bynum '+sel
        currSel=mdUniverse.select_atoms(strSel)

        for ts in mdUniverse.trajectory:
            # particle selection in mdUniverse
            currSel=mdUniverse.select_atoms(strSel)

            # center of mass
            [x,y,z]=currSel.center_of_mass()

            # condition
            condVol = eval(volume)
        
            # add on-off into out array
            if ( condVol == 1) :
                outSel.append(sel)
                outSelLine=np.hstack( (outSelLine , countLine) )
                break

        countLine+=1

    outTuple=(outSelLine,outSel)
        
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

    # particle selection
    strSel='( bynum ' + selection + ' )'

    # out array for bool values
    OnOffArray=np.ones((nFrames), dtype=int)
    
    # iterate over trajectory
    counter=0;
    for ts in mdUniverse.trajectory:
        # particle selection in mdUniverse
        currSel=mdUniverse.select_atoms(strSel)

        # center of mass
        [x,y,z]=currSel.center_of_mass()

        # condition
        condVol = eval(volume)
        
        # add on-off into out array
        if ( condVol == 0) :
            OnOffArray[counter]=-1

        # clean
        currSel=[]
        counter=counter+1

    # out
    return OnOffArray



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


