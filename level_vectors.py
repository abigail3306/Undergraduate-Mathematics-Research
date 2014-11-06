from math import floor
from generate_vectors import generate_level_vectors

# TODO: is there a way to generate vectors more efficently?
# function returns the gap vector of a given level_vector.  It iterates
# through the length of the level vector.
def get_gap_vector(level_vector):
    dtable = []
    for j in range(0,len(level_vector)):
	new_col = []
	for i in range(0,j):
	    entry = level_vector[j] - level_vector[i]
	    if entry < 0:
	       entry = -1 - entry
	       dtable[i].append(entry)
	    else:
	       new_col.append(entry)
	dtable.append(new_col)
    sums = [sum(new_col) for new_col in dtable]
    sums.sort()
    sums.reverse()
    #print 'For level vector', level_vector, ', the gap vector is:',sums
    return [sums[j] - sums[j+1] for j in range(0,len(sums) - 1)]


# function determines the coxeter_length by creating an empty array
# and adding 1 to the position's value if a permutation is found that
# avoids at a particular length.
def coxeter_length(gap_vector):
    offsets = range(1, len(gap_vector) + 1)
    dot_product = [gap_vector[i] * offsets[i] for i in range(0, len(gap_vector))]    
    return sum(dot_product)

# function returns 1-line notation of a given level vector by iterating through
# the values of the level vector and multiplying them by the number of 
# elements + the current row position value
def one_line(level_vector):
    n = len(level_vector)
    ret = []
    for i in range(1,n+1):
	ret.append(level_vector[i-1] * n+i)  
    ret.sort()
    #print 'For level vector', level_vector, ', the one-line notation is:', ret
    return ret


# function determines whether a particular point, i, in the abacus
# is a gap or not and returns a boolean value
def is_gap(level_vector,i):
    n = len(level_vector)
    level_of_i = floor((i-1)/n) 
    level_def_bead = level_vector[(i % n) - 1]

    print 'For level vector', level_vector, ', value', i, 'is', (level_of_i > level_def_bead)
    return(level_of_i > level_def_bead)

