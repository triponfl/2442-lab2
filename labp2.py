#! /bin/env python
# -*- coding: utf-8 -*-

import numpy as np
# Replace these names with your names
authors = ['A. Flavia Tripon', 'B. Agustina Perini Álvarez']

#j is across!!
#i is down!!

# Initiating dynamic programming matrices, S and trace,
# Input is the lengths of each of the two sequences
# Output is the initiated dynamic programing and trace matrices
def initiate_global_dp(m,n):
    S = np.zeros((m+1, n+1))       # An (m+1)*(n+1) matrix, initiated with 0's
    # Careful! S and trace are shape (m+1)*(n+1). The sequences lengths are m and n. Remember that in the future!
    # For the trace matrix we use a three dimentional matrix of booleans
    # where 
    #   trace(x,y,0) indicates a match in x,y
    #   trace(x,y,1) indicates an insert in x,y (fix column)
    #   trace(x,y,2) indicates a delete in x,y (fix row) # so erm actually the 2 is index related. so they are saying there is a T in tuple[2]
    trace = np.zeros((m+1, n+1, 3), dtype=np.bool) # An (m+1)*(n+1)*3 boolean matrix, initiated with (False,False,False)
    #also oops apparently something odd was happening with bool8?? switched to bool
    # why is trace 3 dimensions - scoring you only need to consider row and column and you output a value into that cell. to get the output for trace, you need to consider row, column, and what combination of of T/F (match, insert, gap) you have in that cell. that is the third dimension.
    # 3 is the SIZE of the thrrd axis, not the dimension
    # so (m+1, n+1, 5) would mean that there is a tuple of 5 elements in each cell but still has three dimensions
    

    # HINT: try printing S and trace for debugging.
    # trace is difficult to 3D-visualize out of the box with numpy, so use the function pretty_trace(trace) for a better result!
    # you can also use the function pretty_trace_arrows(trace) for an even clearer (but less abstract) result


    # First initiate the origin of S (0,0) here:
    S[0,0] = 0. # so this sets up the first cell of the scoring matrix - period is to make sure its a float even tho just 0 would work too i think
    trace[0,0,:] = (0.,0.,0.) # this sets up the first cell of the trace matrix. #also the : makes sure it takes every element from cell 0,0. so it calls the whole index of that axis (three in this case)


    # Now, fill in the first row and the first column of the matrices S and trace
    # Initiate the first column of S and trace here:
    for i in range(1,m+1): # +1 b/c we need the empty row/column for the first row/column of the matrix and we did it earlier ?? i think
        S[i,0] = i * gap_penalty() #we use the function below which is already - 2.0 # btw the 0 is not changing because we are initializing. so this is column 0 and we are filling in the cells below ONLY here. not the whole matrix. so i is the reference seq and going down, so if we say there is a gap here, the other sequence MUST have an insert
        trace[i,0,:] = (0.,0.,1) # FFT = true at position 2 (index thing) means gap. needleman is intiialized with gaps so yeah
    
    # Initiate the first row of S and trace here:
    for j in range(1,n+1):
        S[0,j] = j * gap_penalty()
        trace[0,j,:] = (0.,1.,0) # FTF. so because these are different sequences, i thinkwe're saying that while there is a gap in the other one, you can consider this as an insertion, otherwise you're just creating gaps and it'll look like = when you write it out ? 
    # Return the initiated matrices
    return S, trace

# Fill in the dynamic programming matrix and the trace
def global_align(seqA,seqB):
    # Initiating variables
    m, n = len(seqA), len(seqB)  # Be careful with indexing problems, such as +-1 errors. When in doubt, print.
    S,trace = initiate_global_dp(m,n)
    # Fill in the rest of the dynamic programming matrix, and the trace
    for i in range(1,m+1):
        for j in range(1,n+1): #lets up go through the whole matrix. so when we are on row 1 we go through every column (left to right and then down and repeat)
            match_mismatch = S[i-1,j-1] + match_score(seqA[i-1],seqB[j-1])
            insert = S[i,j-1] + gap_penalty() 
            gap = S[i-1,j] + gap_penalty() 
            S[i,j] = max(match_mismatch, insert, gap) #take highest value and make that the score
            #so now to set up the values for following the optimal path
            if match_mismatch >= max(insert, gap):
                trace[i,j,:] = (1.,0.,0.) # TFF for match/mismatch
            elif insert >= max(match_mismatch, gap):
                trace[i,j,:] = (0.,1.,0.) # FTF for insert
            else: # this is gap
                trace[i,j,:] = (0.,0.,1.) # FFT
    score_of_the_alignment = S[m,n] #should be the last cell for needleman
    return S, trace, score_of_the_alignment

# EXCERCISE 1
def initiate_dp(m,n, alignment_type):
    S = np.zeros((m+1, n+1))    
    trace = np.zeros((m+1, n+1, 3), dtype=np.bool) # An (m+1)*(n+1)*3 boolean matrix, initiated with (False,False,False)

    S[0,0] = 0. # so this sets up the first cell of the scoring matrix - period is to make sure its a float even tho just 0 would work too i think
    trace[0,0,:] = (0.,0.,0.) # this sets up the first cell of the trace matrix. #also the : makes sure it takes every element from cell 0,0. so it calls the whole index of that axis (three in this case)

    if alignment_type == 1: #global
        for i in range(1,m+1): # +1 b/c we need the empty row/column for the first row/column of the matrix and we did it earlier ?? i think
            S[i,0] = i * gap_penalty() #we use the function below which is already - 2.0 # btw the 0 is not changing because we are initializing. so this is column 0 and we are filling in the cells below ONLY here. not the whole matrix. so i is the reference seq and going down, so if we say there is a gap here, the other sequence MUST have an insert
            trace[i,0,:] = (0.,0.,1) # FFT = true at position 2 (index thing) means gap. needleman is intiialized with gaps so yeah
        
        for j in range(1,n+1):
            S[0,j] = j * gap_penalty()
            trace[0,j,:] = (0.,1.,0) # FTF. so because these are different sequences, i thinkwe're saying that while there is a gap in the other one, you can consider this as an insertion, otherwise you're just creating gaps and it'll look like = when you write it out ? 
    else:
        for i in range(1, m+1): 
            S[i,0] = 0.
            trace[i,0,:] = (0.,0.,0.)
        for j in range(1, n+1):
            S[0,j] = 0.
            trace[0,j,:] = (0.,0.,0.)
    # Return the initiated matrices
    return S, trace

# Fill in the dynamic programming matrix and the trace
def align(seqA,seqB, alignment_type):
    # Initiating variables
    m, n = len(seqA), len(seqB)  # Be careful with indexing problems, such as +-1 errors. When in doubt, print.
    S,trace = initiate_dp(m,n,alignment_type)
    # Fill in the rest of the dynamic programming matrix, and the trace
    for i in range(1,m+1):
        for j in range(1,n+1): #lets up go through the whole matrix. so when we are on row 1 we go through every column (left to right and then down and repeat)
            match_mismatch = S[i-1,j-1] + match_score(seqA[i-1],seqB[j-1])
            insert = S[i,j-1] + gap_penalty() 
            gap = S[i-1,j] + gap_penalty() 
            if alignment_type == 1:
                S[i,j] = max(match_mismatch, insert, gap) #take highest value and make that the score
            else:
                S[i,j] = max(0, match_mismatch, insert, gap) #local can go back to 0
            #so now to set up the values for following the optimal path
            if alignment_type == 0 and 0>= max(match_mismatch, insert, gap): 
                trace[i,j,:] = (0.,0.,0.) #we're saying that the optimal local alignment restarts at this spot because we don't want any negative values
            elif match_mismatch >= max(insert, gap):
                trace[i,j,:] = (1.,0.,0.) # TFF for match/mismatch
            elif insert >= max(match_mismatch, gap):
                trace[i,j,:] = (0.,1.,0.) # FTF for insert
            else: # this is gap
                trace[i,j,:] = (0.,0.,1.) # FFT
    if alignment_type == 1: #global
        score_of_the_alignment = S[m,n] #should be the last cell for needleman
    else: # local
        score_of_the_alignment = np.max(S)
    return S, trace, score_of_the_alignment #later on the formatting fucniton i think would have to have a specified start_from for local??


#EXCERCISE 2 - really dont know if this is correct and don't know how to check (same with above one tbh)
#gap penalty would be different now
new_gap = -3
cont_gap = -1
#ok so now i need three matrices based off those slides given
# matrix 1 is if our max score ends in match/mismatch
# matrix 2 is if our max score ends with a gap in i
# matric 3 is if our max score ends with a gap in j
# infinity starts coming into play (silde 12) so ig i gotta keep that in mind
#i don't understand
#https://github.com/biopython/biopython/blob/master/Bio/pairwise2.py didn't help that much actually

# Conventions used below (same as the trace layers in the exercises above):
#   Mij[i,j] : best score where a_i is aligned to b_j            (diagonal move)
#   Xij[i,j] : best score where a_i is aligned to a gap          (vertical move, "delete",  trace index 2)
#   Yij[i,j] : best score where b_j is aligned to a gap          (horizontal move, "insert", trace index 1)
# In the affine case each trace layer records WHICH MATRIX the predecessor came from:
#   (T,F,F) = came from M,  (F,T,F) = came from Y,  (F,F,T) = came from X
def initiate_affine_dp(m,n,alignment_type):
    Mij = np.zeros((m+1, n+1)) #not sure if these are correct
    Xij = np.zeros((m+1, n+1))
    Yij = np.zeros((m+1, n+1))
    trace_Mij = np.zeros((m+1, n+1, 3), dtype=np.bool)
    trace_Xij = np.zeros((m+1, n+1, 3), dtype=np.bool)
    trace_Yij = np.zeros((m+1, n+1, 3), dtype=np.bool)

    Mij[0,0] = 0
    # these are other layers so idk if you can initilaize with 0 since this is if there is a gap in seq A or B??
    # FIX: "impossible" is -inf, not None (None becomes NaN in a float array, and max()
    #      with NaN gives silent nonsense). At (0,0) nothing has been consumed, so you
    #      cannot be in the middle of a gap.
    Xij[0,0] = -np.inf
    Yij[0,0] = -np.inf

    if alignment_type == 1: #global
        for i in range(1,m+1):
            # FIX: new_gap and cont_gap are numbers, not functions -> no ()
            Xij[i,0] = new_gap + ((i-1) * cont_gap) #i-1 b/c we need no continued gap at first, just the new gap
            # FIX: the first gap comes from M[0,0]; the following ones extend X
            trace_Xij[i,0,:] = (1,0,0) if i == 1 else (0,0,1)
            # FIX: the rest of column 0 can NOT be reached: you cannot have aligned a_i to
            #      b_0 (M) or be in a horizontal gap (Y) when nothing of B is consumed yet.
            #      Leaving these at 0 lets the recursion "start for free" from the border.
            Mij[i,0] = -np.inf
            Yij[i,0] = -np.inf
        for j in range(1,n+1):
            Yij[0,j] = new_gap + ((j-1) * cont_gap)
            trace_Yij[0,j,:] = (1,0,0) if j == 1 else (0,1,0)
            Mij[0,j] = -np.inf   # FIX: same for row 0
            Xij[0,j] = -np.inf
    else: #local
        for i in range(1, m+1):
            Mij[i,0] = 0.
            trace_Mij[i,0,:] = (0.,0.,0.)
        for j in range(1, n+1):
            Mij[0,j] = 0.
            trace_Mij[0,j,:] = (0.,0.,0.)
        # (Xij/Yij borders can stay 0 here: a local alignment starting with a gap is
        #  never better than starting at the 0 in Mij, so it never wins.)
    # Return the initiated matrices
    return Mij, Xij, Yij, trace_Mij, trace_Xij, trace_Yij

# FIX: small helper so the "which of the three did we come from" logic is written once.
# candidates are given in the order (from M, from Y, from X) so that the argmax
# index is directly the trace layer index.
def _best(from_M, from_Y, from_X):
    cands = (from_M, from_Y, from_X)
    k = int(np.argmax(cands))
    tr = [0,0,0]
    tr[k] = 1
    return cands[k], tuple(tr)

def affine_align(seqA,seqB, alignment_type):
    # Initiating variables
    m, n = len(seqA), len(seqB)
    Mij, Xij, Yij,trace_Mij, trace_Xij, trace_Yij = initiate_affine_dp(m,n,alignment_type)
    for i in range(1,m+1):
        for j in range(1,n+1):
            # M: a_i aligned to b_j. Predecessor is the diagonal cell in ANY of the three matrices.
            # FIX: match_score is a function -> () not []
            s = match_score(seqA[i-1],seqB[j-1])
            Mij[i,j], tr = _best(Mij[i-1,j-1] + s, Yij[i-1,j-1] + s, Xij[i-1,j-1] + s)
            if alignment_type == 0 and Mij[i,j] <= 0: #local can restart here
                Mij[i,j] = 0.
                tr = (0,0,0)
            # FIX: assign INTO the array (trace_Mij[i,j,:] = ...), not trace_Mij = ...
            #      (the old code replaced the whole matrix by a tuple)
            trace_Mij[i,j,:] = tr

            # X: a_i aligned to a gap.
            # FIX: the predecessor is the cell ABOVE, [i-1, j] -- not [i-1, j-1].
            #      A vertical gap only consumes a character from seqA.
            Xij[i,j], trace_Xij[i,j,:] = _best(Mij[i-1,j] + new_gap,
                                              Yij[i-1,j] + new_gap,
                                              Xij[i-1,j] + cont_gap)

            # Y: b_j aligned to a gap.
            # FIX: the predecessor is the cell to the LEFT, [i, j-1].
            # FIX: the old global branch filled Yij with X's candidates (copy/paste).
            Yij[i,j], trace_Yij[i,j,:] = _best(Mij[i,j-1] + new_gap,
                                              Yij[i,j-1] + cont_gap,
                                              Xij[i,j-1] + new_gap)

    if alignment_type == 1: #global
        score_of_the_alignment = max(Mij[m,n],Xij[m,n],Yij[m,n]) #should be the last cell for needleman
    else: # local
        # FIX: np.max(a, b, c) means np.max(a, axis=b, out=c). Local score = best cell in M.
        score_of_the_alignment = np.max(Mij)
    # FIX: the score was computed but never returned
    return Mij, Xij, Yij, trace_Mij, trace_Xij, trace_Yij, score_of_the_alignment

# FIX (new): traceback for the affine case. format_alignment() above only knows one trace
# layer; here we also have to remember which matrix we are currently in.
def format_affine_alignment(seqA, seqB, M, X, Y, trace_M, trace_X, trace_Y, start_from=None):
    if start_from:                 # local: (i, j) of the best cell in M
        i, j = start_from
        cur = 0
    else:                          # global: end cell, in whichever matrix holds the best score
        i, j = len(seqA), len(seqB)
        cur = int(np.argmax((M[i,j], Y[i,j], X[i,j])))   # 0 = M, 1 = Y, 2 = X
    traces = (trace_M, trace_Y, trace_X)
    outA, outB = "", ""
    while i > 0 or j > 0:
        tr = traces[cur][i,j]
        if not tr.any():        # local alignment start
            break
        nxt = int(np.argmax(tr))   # which matrix the predecessor is in
        if cur == 0:               # in M: consume both
            i, j = i-1, j-1
            outA = seqA[i] + outA
            outB = seqB[j] + outB
        elif cur == 1:             # in Y: gap in A, consume b_j
            j = j-1
            outA = "-" + outA
            outB = seqB[j] + outB
        else:                      # in X: gap in B, consume a_i
            i = i-1
            outA = seqA[i] + outA
            outB = "-" + outB
        cur = nxt
    return outA, outB

#sources
#the textbook given and that's about it
#https://www.youtube.com/watch?v=DQQ_q2dn2ds
#https://www.youtube.com/watch?v=um8h3P216Fk
#https://github.com/biopython/biopython/blob/master/Bio/pairwise2.py

## The following functioins are give to you as a help
# Return the gap penalty
def gap_penalty():
    return -2.0

# Return the match score of letterA and letterB.
# If one of the letters is a gap, return the gap penalty
# otherwise return their match/mismatch score
def match_score(letterA,letterB):
    if letterA == '-' or letterB == '-':
        return gap_penalty()
    elif letterA == letterB:
        return 3.0
    else:
        return -1.0
    
# Print 2 sequences on top of each other
def print_alignment(seqA,seqB):
    print(seqA)
    print(seqB)

# Print a dynamic programming score matrix
# together with its sequences
def print_dynamic(seqA,seqB,dpm):
    seqA,seqB = "-" + seqA, "-" + seqB
    m,n = len(seqA),len(seqB)
    print('{:^5}'.format(' '), end=''),
    for j in range(n):
        print('{:^5}'.format(seqB[j]), end='')
    print()
    for i in range(m):
        print ('{:^5}'.format(seqA[i]), end="")
        for j in range(n):
            print ('{:5.1f}'.format(dpm[i,j]), end="")
        print()
    print()

def pretty_trace(trace):
    m, n, _ = trace.shape
    str_cells = [[str(tuple(trace[i, j])).replace("True", "T").replace("False", "F")
                  for j in range(n)] for i in range(m)]
    col_widths = [max(len(str_cells[i][j]) for i in range(m)) for j in range(n)]
    for i in range(m):
        row_str = " ".join(f"{str_cells[i][j]:>{col_widths[j]}}" for j in range(n))
        print(row_str)

def pretty_trace_arrows(trace):
    m, n, _ = trace.shape
    arrows = {0: "↖", 1: "←", 2: "↑"}
    str_cells = []
    for i in range(m):
        row = []
        for j in range(n):
            dirs = "".join(arrows[k] for k in range(3) if trace[i, j, k])
            row.append(dirs if dirs else ".")
        str_cells.append(row)
    col_widths = [max(len(str_cells[i][j]) for i in range(m)) for j in range(n)]
    total_width = sum(col_widths) + (n - 1) * 3
    lines = ["-" * total_width]
    for i in range(m):
        row_str = " | ".join(f"{str_cells[i][j]:>{col_widths[j]}}" for j in range(n))
        lines.append(row_str)
    lines.append("-" * total_width)
    return "\n".join(lines)


# Format an alignment by inserting gaps in sequences given a trace matrix
def format_alignment(seqA, seqB, trace, start_from = None):
    if start_from:
        i, j = start_from
    else:
        i, j = len(seqA), len(seqB)
    outA, outB = "",""
    while i>0 or j>0:
        if trace[i,j,0]: # match
            i, j = i-1, j-1
            outA = seqA[i] + outA
            outB = seqB[j] + outB
        elif trace[i,j,1]: # insert
            i, j = i, j-1
            outA = "-" + outA
            outB = seqB[j] + outB
        elif trace[i,j,2]: # delete
            i, j = i-1, j
            outA = seqA[i] + outA
            outB = "-" + outB
    return outA,outB



# Test code for the dna2aa function. 
# Will only be executed if this file is run directly
# e.g. by running the command "python labp2.py"
if __name__ == "__main__":
    seqA, seqB = "ATG", "GAT"
    dp, trace, max_score = global_align(seqA, seqB)
    print_dynamic(seqA, seqB, dp)
    pretty_trace(trace) 
    pretty_trace_arrows(trace) 
    #print('\n'.join(format_alignment(seqA, seqB, trace)))
    print(f"Score: {max_score}")

    # FIX (new): a small test of the affine version. With open=-3, extend=-1 one long
    # gap should be preferred over several short ones.
    seqA, seqB = "GATTACA", "GCA"
    M, X, Y, tM, tX, tY, score = affine_align(seqA, seqB, 1)
    print_dynamic(seqA, seqB, M)
    print('\n'.join(format_affine_alignment(seqA, seqB, M, X, Y, tM, tX, tY)))
    print(f"Affine global score: {score}")   # expected 3.0: G,C,A matched (+9), one gap of 4 (-3-1-1-1)
