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
def initiate_affine_dp(m,n,alignment_type):
    Mij = np.zeros((m+1, n+1)) #not sure if these are correct
    Xij = np.zeros((m+1, n+1)) 
    Yij = np.zeros((m+1, n+1))
    trace_Mij = np.zeros((m+1, n+1, 3), dtype=np.bool)
    trace_Xij = np.zeros((m+1, n+1, 3), dtype=np.bool)
    trace_Yij = np.zeros((m+1, n+1, 3), dtype=np.bool)

    Mij[0,0] = 0
    # these are other layers so idk if you can initilaize with 0 since this is if there is a gap in seq A or B??
    Xij[0,0] = None
    Yij[0,0] = None

    if alignment_type == 1: #global
        for i in range(1,m+1): 
            Xij[i,0] = new_gap() + ((i-1) * cont_gap()) #i-1 b/c we need no continued gap at first, just the new gap 
            trace_Xij[i,0,:] = (0.,0.,1) # FFT = true at position 2 (index thing) means gap. needleman is intiialized with gaps so yeah
        for j in range(1,n+1):
            Yij[0,j] = new_gap() + ((j-1) * cont_gap())
            trace_Yij[0,j,:] = (0.,1.,0) #just the same as above ig? imagining new gap and then being continued just horizontally? (or insertions?)
            #ot sure where Mij comes into play here??
    else: #local
        for i in range(1, m+1): 
            Mij[i,0] = 0.
            trace_Mij[i,0,:] = (0.,0.,0.)
        for j in range(1, n+1):
            Mij[0,j] = 0.
            trace_Mij[0,j,:] = (0.,0.,0.)
    # Return the initiated matrices
    return Mij, Xij, Yij, trace_Mij, trace_Xij, trace_Yij
# apparently needle is the only one who initializes corrrectly and i cant figure out how - embAlignPathCalcWithEndGapPenalties

def affine_align(seqA,seqB, alignment_type):
    # Initiating variables
    m, n = len(seqA), len(seqB) 
    Mij, Xij, Yij,trace_Mij, trace_Xij, trace_Yij = initiate_affine_dp(m,n,alignment_type)
    for i in range(1,m+1):
        for j in range(1,n+1): #lets up go through the whole matrix. so when we are on row 1 we go through every column (left to right and then down and repeat)
            #three matrices/routes to consider?

            #so if previous was a match/mismatch - M represents match/mismatch so it is the focus here?
            previous = (Mij[i-1,j-1], Xij[i-1,j-1], Yij[i-1,j-1]) #calling all previous possibilities
            best_previous = max(previous)
            Mij[i,j] = best_previous + match_score[seqA[i-1],seqB[j-1]]
            
            if alignment_type == 0: #local
                Mij[i,j] = max(0., best_previous + match_score[seqA[i-1],seqB[j-1]])
            elif alignment_type == 0 and Mij[i,j] <= 0: 
                trace_Mij[i,j,:] = (0.,0.,0.)
            elif alignment_type == 1 and previous[0] >= max(previous[1], previous[2]): # global - if new match/mismatch
                trace_Mij = (1.,0.,0.)
            elif alignment_type == 1 and previous[1] >= max(previous[0], previous[2]): # global - if new gap
                trace_Mij = (0.,0.,1.)
            elif alignment_type == 1 and previous[2] >= max(previous[0], previous[1]): # global - if new insert
                trace_Mij = (0.,1.,0.)
            else:
                print('something went wrong M')

            #Xij now - previous was a gap (in seqA -> i)
            Xij_ext_A = Xij[i-1,j-1] + cont_gap
            #others would be new gap
            Mij_new_A = Mij[i-1,j-1] + new_gap
            Yij_new_A = Yij[i-1,j-1] + new_gap
            #so we have the three possible routes now?
            Xij[i,j] = max(Xij_ext_A, Mij_new_A, Yij_new_A)

            if alignment_type == 0: #local
                Xij[i,j] = max(0., max(Xij_ext_A, Mij_new_A, Yij_new_A))
            elif alignment_type == 0 and Xij[i,j] <= 0: 
                trace_Xij[i,j,:] = (0.,0.,0.)
            elif alignment_type == 1 and Xij_ext_A >= max(Mij_new_A, Yij_new_A): # global - if new match/mismatch
                trace_Xij = (0.,1.,0.)
            elif alignment_type == 1 and Mij_new_A >= max(Xij_ext_A, Yij_new_A): # global - if new gap
                trace_Xij = (1.,0.,0.)
            elif alignment_type == 1 and Yij_new_A >= max(Mij_new_A, Xij_ext_A): # global - if new insert
                trace_Xij = (0.,0.,1.)
            else:
                print('something went wrong X')

            #Yij now - previous was an insert (in seqA -> i) - so gap in B
            Yij_ext_B = Yij[i-1,j-1] + cont_gap
            #others would be new gap
            Mij_new_B = Mij[i-1,j-1] + new_gap
            Xij_new_B = Xij[i-1,j-1] + new_gap
            #so we have the three possible routes now?
            Yij[i,j] = max(Xij_ext_A, Mij_new_A, Yij_new_A)

            if alignment_type == 0: #local
                Yij[i,j] = max(0., max(Yij_ext_B, Mij_new_B, Xij_new_B))
            elif alignment_type == 0 and Xij[i,j] <= 0: 
                trace_Yij[i,j,:] = (0.,0.,0.)
            elif alignment_type == 1 and Yij_ext_B >= max(Mij_new_B, Xij_new_B): # global - if new match/mismatch
                trace_Yij = (0.,1.,0.)
            elif alignment_type == 1 and Mij_new_B >= max(Yij_ext_B, Xij_new_B): # global - if new gap
                trace_Yij = (1.,0.,0.)
            elif alignment_type == 1 and Xij_new_B >= max(Mij_new_B, Yij_ext_B): # global - if new insert
                trace_Yij = (0.,0.,1.)
            else:
                print('something went wrong Y')

    if alignment_type == 1: #global
        score_of_the_alignment = max(Mij[m,n],Xij[m,n],Yij[m,n]) #should be the last cell for needleman
    else: # local
        score_of_the_alignment = np.max(Mij[m,n],Xij[m,n],Yij[m,n])
    return Mij, Xij, Yij, trace_Mij, trace_Xij, trace_Yij #no idea how to test this
 

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
