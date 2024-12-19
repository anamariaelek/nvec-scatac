import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import pickle as pk
from tqdm import tqdm
from numba import jit, njit, prange

# Function to calculate GC content of a fasta file
# @param path: The path to a fasta-formatted genome
# @return fraction of all bases that are G or C
def calculateGenomeGCContent( path ):
    gc = 0
    at = 0
    for record in SeqIO.parse( path, 'fasta' ):
        s = record.seq.upper()
        gc = gc + s.count('G')
        gc = gc + s.count('C')
        at = at + s.count('A')
        at = at + s.count('T')
    return gc / (gc+at)

# Function to read in MoDISco PFMs from the output directory
# @param path: the base directory to search for pfms
#   should contain MoDISco outputs from interpreting a single model
# @param task: which task (counts or profile) interpretation
#   to read PFMs from
# @param prefix: prefix to add to the PFM names
def readPFMsFromFolder( path, task='profile', prefix=None ):
    task = task.lower()
    if task not in [ 'profile', 'counts' ]:
        print( 'Invalid task: choose either \'profile\' or \'counts\'' )
        return None
    p = prefix
    if prefix is None:
        p = str(path).strip('/').split('/')[-1]
    pfmpath = os.path.join( path, 'pfms_{0}.npy'.format(task) )
    countpath = os.path.join( path, 'pfms_{0}_n_seqlets.npy'.format(task) )
    inpfms = np.load( pfmpath, allow_pickle=True ).flat[0]
    incounts = np.load( countpath, allow_pickle=True ).flat[0]
    names = []
    outcounts = []
    outpfms = {}
    for pattern in inpfms:
        name = '_'.join([p,task,pattern])
        names.append(name)
        outcounts.append(incounts.get(pattern))
        outpfms[name] = inpfms[pattern]
    return outpfms, names, outcounts

# Generate a one-hot encoded PFM
# from a string of A,T,G,C
# Output will be a N x 4 numpy array
def pfmFromSeq( seq ):
    s = seq.upper()
    baseMap = { 'A':0, 'C':1, 'G':2, 'T':3 }
    pfm = np.zeros((len(seq),4))
    for i in range(len(seq)):
        base = s[i]
        pos = baseMap.get(base)
        if pos is None:
            return None
        pfm[i,pos] = 1
    return pfm

# @param pfm: an (N,4) numpy array representing a pfm
# @param pseudo: floating point value representing
#   the pseudocount to add to the pfm (avoids log(0))
# @return the pfm with pseudocount added to all entries
#   and then normalized such that rows still sum to 1
def addPFMPseudoCount( pfm, pseudo ):
    plusPseudo = pfm + pseudo
    return plusPseudo / plusPseudo.sum(1).reshape(-1,1)

# @param pfm: an (N,4) numpy array representing a pfm
# @param bkGC: GC content for the genome of origin
# @param pseudo: pseudocount to add to the pfm
# @return pfm weighted by the information content of 
#   each base relative to the background GC content
#   frequencies above background will be positive,
#   frequencies below background will be negative
def relativeInfo( pfm, bkgGC=0.5, pseudo=0 ):
    gc = bkgGC / 2
    at = 0.5 - gc
    background = np.array([[at,gc,gc,at]])
    ppfm = addPFMPseudoCount( pfm, pseudo )
    return ppfm*(np.log2(ppfm)-np.log2(background)) + 0

# @param pfm: an (N,4) numpy array representing a pfm
# @param cutoff: minimum total information content per
#   position for trimming the PFM
# @param bkGC: GC content for the genome of origin
# @return the pfm trimmed between the first and last
#   positions with total information content above the cutoff
def trimPFM( pfm, cutoff, bkgGC=0.5 ):
    info = relativeInfo( pfm, bkgGC ).sum(1)
    goodInd = np.argwhere( info >= cutoff ).flatten()
    trimmed = pfm.copy()
    if len(goodInd) > 0:
        start = goodInd.min()
        end = goodInd.max() + 1
        trimmed = trimmed[start:end,:]
    return trimmed

# Generate a logomaker Logo object from a PFM
# @param pfm: an (N,4) numpy array representing a pfm
# @param pseudo: floating point value representing
#   the pseudocount to add to the pfm (avoids log(0))
def makeLogoFromPFM( pfm, pseudo=0.001, bkgGC=0.5 ):
    info = relativeInfo( addPFMPseudoCount(pfm,pseudo), bkgGC )
    df = pd.DataFrame( data=info, columns=['A','C','G','T'] )
    logo = logomaker.Logo( df )
    logo.ax.set_yticks([0,1,2])
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    return logo

# @param A,B: a pair of (N,4) numpy arrays 
#   representing two pfms to compare
# @return the result of taking the Jensen-Shannon
#   distance between every pair of rows, output
#   as a length-N numpy array
def pairwiseJSD( A, B ):
    if A.shape != B.shape:
        print( 'Shapes of input arrays do not match!' )
        return np.array([-1.0])
    M = (A+B)/2
    JSD = A*(np.log2(A)-np.log2(M)) + \
          B*(np.log2(B)-np.log2(M))
    # sum it up and aggregate
    JSD = JSD.sum(1) / 2
    return JSD

# Get the min-JSD ungapped alignment of two PFMs
# A psuedocount of pseudo will be added to all bases
# to avoid taking log(0)
# The total JSD will be aggregated over ...
#   - All entries of the input pair for overlap='complete' (default)
#   - All entries of the first PFM for overlap='fixed'
#   - All entries of the shorter PFM for overlap='min'
# JSD values in this region can aggregate through ...
#   - Summation for agg='sum' (default)
#   - Averaging for agg='average'
#   - Maximum value for agg='max'
#   - Softmax-weighted average based on information content
#     for agg='infoAverage'
# Per-base JSD values below jsdCutoff will be pushed to 0
# Only alignments with at least minOverlap aligning bases
# will be considered in the output. In the case that either
# input PWM is smaller than minOverlap, the total length of
# the two input PWMs will be output instead
# @param fixed: (N,4) numpy array representing the pfm  
#   that will be the basis for comparison
# @param shifted: (N,4) numpy array representing the pfm 
#   that will be shifted and rotated to compare to fixed
# @param overlap: which positions to aggregate JSD over (described above)
# @param pseudo: pseudocount to add to all positions
# @param agg: method of aggregating JSD values (described above)
# @param jsdCutoff: minimum JSD value per-position to consider meaningful
# @param minOverlap: minimum number of overlapping positions to compare
# @param bkGCF: background GC content for fixed pfm
# @param bkGCS: background GC content for shifted pfm
# @param removeBackground: whether to subtract differences in base
#   frequencies due to background GC content in the two PFMs
# @return 3-tuple of ( best aggregate JSD value over all possible shifts,
#                      starting position of shifted pfm relative to fixed
#                      pfm for the best possible alignment,
#                      whether or not the best alignment requires the 
#                     reverse complement of the shifted PFM )
def bestOverlapPFMs( fixed, shifted, overlap='complete', pseudo=1e-3, 
                     agg='sum', jsdCutoff=0, minOverlap=4,
                     bkgGCF=0.5, bkgGCS=None, removeBackground=False ):
    # check arguments are valid
    if overlap not in [ 'complete', 'fixed', 'min' ]:
        overlap = 'complete'
    if agg not in [ 'sum', 'average', 'max', 'infoAverage' ]:
        agg = 'sum'
    if bkgGCS is None:
        bkgGCS = bkgGCF
        
    f = addPFMPseudoCount( fixed, pseudo )
    s = addPFMPseudoCount( shifted, pseudo )
    
    # check that both pfms are above the minimum overlap
    lenF = len(f)
    lenS = len(s)
    if lenS < minOverlap or lenF < minOverlap:
        return (lenS + lenF,0,False)
    
    gcF = bkgGCF / 2
    atF = 0.5 - gcF
    backgroundF = np.array([[atF,gcF,gcF,atF]])
    # array for storing scores at each overlap
    overlaps = np.zeros(lenF+lenS+1-2*minOverlap)
    # extended version of fixed, padded with background
    extF = np.ones((2*lenS+lenF-2*minOverlap,4))
    extF[:,0] = extF[:,3] = atF
    extF[:,1] = extF[:,2] = gcF
    extF[lenS-minOverlap:lenF+lenS-minOverlap,:] = f
    gcS = bkgGCS / 2
    atS = 0.5 - gcS
    backgroundS = np.array([[atS,gcS,gcS,atS]])
    # compute difference in background frequencies 
    # if we want to remove the background
    M = (backgroundF+backgroundS) / 2
    deltaBackground = backgroundF*(np.log2(backgroundF)-np.log2(M)) + \
                      backgroundS*(np.log2(backgroundS)-np.log2(M))
    deltaBackground = deltaBackground.sum() / 2
    for i in prange(lenF+lenS+1-2*minOverlap):
        # shift the PWM relative to the fixed one
        extS = np.ones_like(extF)
        extS[:,0] = extS[:,3] = atS
        extS[:,1] = extS[:,2] = gcS
        extS[i:i+lenS,:] = s
        # find the average of the two distributions
        M = (extF+extS)/2
        # calculate individual terms of the JSD
        JSD = extF*(np.log2(extF)-np.log2(M)) + \
              extS*(np.log2(extS)-np.log2(M))
        # sum it up and aggregate
        JSD = JSD.sum(1) / 2
        if overlap == 'min' and lenF >= lenS:
            start = i
            end = i+lenS
        elif overlap == 'fixed' or overlap == 'min':
            start = lenS-minOverlap
            end = lenS+lenF-minOverlap
        else:
            start = min(i,lenS-minOverlap)
            end = max(lenS-minOverlap+lenF,i+lenS)
        start = int(start)
        end = int(end)
        JSD = JSD[start:end]
        cut = JSD >= jsdCutoff
        if agg == 'sum':
            overlaps[i] = JSD[cut].sum()
            if removeBackground:
                overlaps[i] = overlaps[i] - (end-start) * deltaBackground
        elif agg == 'average':
            overlaps[i] = JSD[cut].sum() / (end-start)
            if removeBackground:
                overlaps[i] = overlaps[i] - deltaBackground
        elif agg == 'infoAverage':
            iF = extF*(np.log2(extF)-np.log2(backgroundF)) + 0
            iS = extS*(np.log2(extS)-np.log2(backgroundS)) + 0
            w = iF.sum(1) + iS.sum(1)
            w = np.exp(w[start:end])[cut]
            w = w / w.sum()
            overlaps[i] = (w*JSD[cut]).sum()
        elif agg == 'max':
            overlaps[i] = JSD.max()
            if removeBackground:
                overlaps[i] = overlaps[i] - deltaBackground
        
    bestOverlap = overlaps.min()
    bestOffset = np.argmin(overlaps) - lenS + minOverlap
    
    # Repeat with the reverse complement
    rc = s[::-1,::-1]
    overlaps = np.zeros(lenF+lenS+1-2*minOverlap)
    for i in prange(lenF+lenS+1-2*minOverlap):
        # shift the PWM relative to the fixed one
        extS = np.ones_like(extF)
        extS[:,0] = extS[:,3] = atS
        extS[:,1] = extS[:,2] = gcS
        extS[i:i+lenS,:] = rc
        M = (extF+extS)/2
        JSD = extF*(np.log2(extF)-np.log2(M)) + \
              extS*(np.log2(extS)-np.log2(M))
        # sum it up and aggregate
        JSD = JSD.sum(1) / 2
        if overlap == 'min' and lenF >= lenS:
            start = i
            end = i+lenS
        elif overlap == 'fixed' or overlap == 'min':
            start = lenS-minOverlap
            end = lenS+lenF-minOverlap
        else:
            start = min(i,lenS-minOverlap)
            end = max(lenS-minOverlap+lenF,i+lenS)
        start = int(start)
        end = int(end)
        JSD = JSD[start:end]
        cut = JSD >= jsdCutoff
        if agg == 'sum':
            overlaps[i] = JSD[cut].sum()
            if removeBackground:
                overlaps[i] = overlaps[i] - (end-start) * deltaBackground
        elif agg == 'average':
            overlaps[i] = JSD[cut].sum() / (end-start)
            if removeBackground:
                overlaps[i] = overlaps[i] - deltaBackground
        elif agg == 'infoAverage':
            iF = extF*(np.log2(extF)-np.log2(backgroundF)) + 0
            iS = extS*(np.log2(extS)-np.log2(backgroundS)) + 0
            w = iF.sum(1) + iS.sum(1)
            w = np.exp(w[start:end])[cut]
            w = w / w.sum()
            overlaps[i] = (w*JSD[cut]).sum()
        elif agg == 'max':
            overlaps[i] = JSD.max()
            if removeBackground:
                overlaps[i] = overlaps[i] - deltaBackground
        
    doRC = False
    if overlaps.min() < bestOverlap:
        doRC = True
        bestOverlap = overlaps.min()
        bestOffset = np.argmin(overlaps) - lenS + minOverlap
        
    return (bestOverlap,bestOffset,doRC)

# @param A: an (N1,4) numpy array representing a pfm
# @param B: an (N2,4) numpy array representing a pfm
# @param shift: number of positions to shift the start 
#   of B relative to the start of A
# @param RCB: whether to take the reverse complement
#   of B when generating the alignment
# @return 4-tuple of ( (N,4) numpy array representing a 
#                      padded version of A,
#                      (N,4) numpy array representing a
#                      padded, shifted, and possibly
#                      reverse-complemented version of B
#                      index of first position where both
#                      pfms are above background,
#                      index of last position where both
#                      pfms are above background
def alignPFMPair( A, B, shift, RCB ):
    ovlStart = np.abs(shift)
    # Non-negative shift -> B starts to the right of A
    if shift >= 0:
        # find width of the alignment
        width = max(A.shape[0],shift+B.shape[0])
        # and end position for the overlap
        ovlEnd = min(A.shape[0],ovlStart+B.shape[0])
        # depends on whether A is longer than B after
        # shifting (first options), or not (second options)
    # Negative shift -> B starts to the left of A
    else:
        # Similar deal to above but reversed
        width = max(B.shape[0],-shift+A.shape[0])
        ovlEnd = min(B.shape[0],ovlStart+A.shape[0])
    if RCB:
        b = B[::-1,::-1]
    else:
        b = B.copy()
        
    # put A and B into padded arrays long 
    # enough for the aligned pfm pairs
    extA = 0.25*np.ones((width,4))
    extB = 0.25*np.ones((width,4))
    if shift >= 0:
        extA[:A.shape[0],:] = A
        extB[shift:shift+B.shape[0],:] = b
    else:
        extB[:B.shape[0],:] = b
        extA[-shift:-shift+A.shape[0]] = A
    return (extA, extB, ovlStart, ovlEnd)

# Generate a one-hot encoded PFM
# from a string of A,T,G,C
# Output will be a N x 4 numpy array
def pfmFromSeq( seq ):
    s = seq.upper()
    baseMap = { 'A':0, 'C':1, 'G':2, 'T':3 }
    pfm = np.zeros((len(seq),4))
    for i in range(len(seq)):
        base = s[i]
        pos = baseMap.get(base)
        if pos is None:
            return None
        pfm[i,pos] = 1
    return pfm

# Align a pair of pfms with know shift and 
#   reverse complementation, then average all positions
# @param A,B,shift,RCB as in alignPFMPair
# @param w1, wb: floating point weights for the 
#   relative contribution of each pfm to the average
# @return weighted average of A and B after shifting
#   and possibly reverse-complementing
def averagePWMs( A, B, shift, RCB, wA=1, wB=1 ):
    # align and shift the pfms as above
    if shift >= 0:
        width = max(A.shape[0],shift+B.shape[0])
    else:
        width = max(B.shape[0],-shift+A.shape[0])
    if RCB:
        b = B[::-1,::-1]
    else:
        b = B.copy()
        
    extA = 0.25*np.ones((width,4))
    extB = 0.25*np.ones((width,4))
    if shift >= 0:
        extA[:A.shape[0],:] = A
        extB[shift:shift+B.shape[0],:] = b
    else:
        extB[:B.shape[0],:] = b
        extA[-shift:-shift+A.shape[0]] = A
    
    # return the weighted average    
    return ( wA*extA + wB*extB )/(wA+wB)

# Function to generate a multiple alignment of several
#   pfms, anchored to a single base pfm
# @param base: an (N1,4) numpy array representing 
#   the pfm to align all others to
# @param others: a list of M (ni,4) numpy arrays
#   representing the pfms to be aligned to base
# @param bkgGC: GC content to pad the arrays with
# @return (N,4,M+1) numpy array containing all of
#   the aligned pfms along the last axis
def generatePFMAlignment( base, others, bkgGC=0.5 ):
    if not isinstance( others, list ):
        print( 'Input must be a list of PFMs!!!!' )
        return base.reshape(-1,4,1)
    # Find best alignments with the base pfm
    starts = np.zeros(len(others),dtype=int)
    ends = np.zeros_like(starts,dtype=int)
    rcs = np.zeros_like(starts,dtype=bool)
    for i in prange(len(others)):
        pfm = others[i]
        _, s, r = bestOverlapPFMs( base, pfm, overlap='complete', 
                                   agg='sum', bkgGCF=bkgGC )
        starts[i] = s
        ends[i] = s + pfm.shape[0]
        rcs[i] = r 
    # Create the alignment of shape N_bases x 4 x N_pfms
    width = int(max(base.shape[0],ends.max()) - min(0,starts.min()))
    baseStart = max(0,-starts.min())
    gc = bkgGC / 2
    at = 0.5 - gc
    align = np.ones((width,4,len(others)+1))
    align[:,0] = align[:,3] = at
    align[:,1] = align[:,2] = gc
    align[baseStart:baseStart+base.shape[0],:,0] = base
    for i in prange(len(others)):
        pfm = others[i]
        if rcs[i]:
            pfm = pfm[::-1,::-1]
        align[baseStart+starts[i]:baseStart+ends[i],:,i+1] = pfm
    return align

# Compute windowed sums over an array
# @param A: a length-N numpy array
# @param windowSize: int size of contiguous windows
#   to sum the array over
# @return a length-(N-windowSize+1) numpy array 
#   containing all of the windowed sums
def windowSum( A, windowSize ):
    csum = np.cumsum( A )
    wsum = csum[windowSize-1:]
    wsum[1:] = wsum[1:] - csum[:-windowSize]
    return wsum

# Find the maximum difference between two pfms over a set window
# @param A: an (N1,4) numpy array representing a pfm
# @param B: an (N2,4) numpy array representing a pfm
# @param windowSize: the integer size of windows to sum over
def getMaxWindowedDistance( A, B, windowSize ):
    # make the alignment
    align = generatePFMAlignment( A, [B] )
    # if the alignment is smaller than the window, pad it
    if align.shape[0] < windowSize:
        align = np.pad( align, ((0,windowSize-align.shape[0]),
                                (0,0),(0,0)),
                        constant_values=0.25 )
    # add pseudocounts to avoid log(0)
    extA = addPFMPseudoCount( align[:,:,0], 0.001 )
    extB = addPFMPseudoCount( align[:,:,1], 0.001 )
    # find JSD per-position
    M = (extA + extB)/2
    JSD = extA*(np.log2(extA)-np.log2(M)) + \
          extB*(np.log2(extB)-np.log2(M))
    JSD = JSD.sum(1) / 2
    # add up the JSD values within all windows and
    # return the max
    wsums = windowSum( JSD, windowSize )
    return wsums.max()

# Function to write a set of pfms to a file in meme format
# @param pfmDict: a dictionary mapping pfm names to 
#   their numpy array representations
# @param fname: the output file name as a string
# @param order: an iterable indicating the order to 
#   write the pfms in (optional)
# @param counts: the number of seqlets supporting
#   each pfm; must be in the same order as order
#   and can only be provided if order is provided
def writeMEME( pfmDict, fname, order=None, counts=None ):
    header = 'MEME version 4\n\n' + \
             'ALPHABET= ACGT\n\n' + \
             'strands: + -\n\n' + \
             'Background letter frequencies\n' + \
             'A 0.25 C 0.25 G 0.25 T 0.25\n\n'
    if order is None:
        order = list(pfmDict.keys())
    if counts is None:
        counts = [20]*len(order)
    with open(fname,'w') as f:
        f.write( header )
        for i, mname in enumerate(order):
            pfm = pfmDict[mname]
            f.write( 'MOTIF {0} {0}\n'.format(mname) )
            f.write( 'letter-probability matrix: alength= 4' +\
                     ' w= {0} nsites= {1} E= 0\n'.format(pfm.shape[0],
                                                         counts[i]) )
            f.write( str(pfm).replace('[','')\
                             .replace(']','')\
                             .replace('\n ','\n')\
                             .replace('0. ','0 ')\
                             .replace('1. ','1 ')\
                             .replace('0.\n','0\n')\
                             .replace('1.\n','1\n')\
                             .strip()\
                             .strip('.') )
            f.write('\n\n')
        f.close()

# Function to read motifs from a meme format file
# @param fname: the path to the meme file as a string
# @return 2-tuple of ( dictionary mapping pfm names to (N,4)
#                      numpy array representations,
#                      list of pfm names in the order they
#                      were written in the file
def readMEME( fname ):
    import numpy as np
    outPFMs = {}
    outNames = []
    pfm = []
    name = ''
    with open( fname, 'r' ) as f:
        for line in f.readlines():
            l = line.strip()
            if len(l) == 0:
                if len(pfm) > 0:
                    outNames.append(name)
                    outPFMs[name] = np.array(pfm).astype('float')
                    pfm = []
                continue
                
            if l.startswith('MOTIF'):
                name = l.split()[1]
            elif l[0].isdigit():
                pfm.append( l.split()[:4] )
            elif len(pfm) > 0:
                outNames.append(name)
                outPFMs[name] = np.array(pfm).astype('float')
                pfm = []
        if len(pfm) > 0:
            outNames.append(name)
            outPFMs[name] = np.array(pfm).astype('float')
            pfm = []
    return outPFMs, np.array(outNames)