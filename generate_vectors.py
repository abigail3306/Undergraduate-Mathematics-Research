from copy import copy
import math
import sys

# iterative + incrementive function that acts as an odometer and gets the next
# value and iterates up through the elements until the max_value of elements
# on that current row are met. This is used to generate gap vectors and level
# vectors
def get_next_row(current_row, value):
    next_row = copy(current_row)
    j = len(next_row) - 2
        
    while (next_row[j] >=  value):
	next_row[j] = -value
	j = j - 1
  
    next_row[j] = next_row[j] + 1
    setLastValue(next_row)
    return next_row


# setLastValue finds the -sum of the row up to the 2nd to last element.  This
# zeros out the row so that the permutation is correct.
def setLastValue(next_row):
    last = -(sum(next_row[0:-1]))
    next_row[-1] = last


# function checks if the row is the lowest row possible. 
# ex. [-cutoff,-cutoff,-cutoff] because this is the row we iniialize with
# the row's last value to equal the sum to 0 does not get set in the original
# function so we set it here.
def setMin(current_row, cutoff):
    if all(current_row[0] == cutoff for cutoff in current_row):	
	row = get_next_row(current_row, cutoff)
        next_row = copy(current_row)

	setLastValue(next_row)


# helper function for get_next_row that adds the rows to a list of lists.  When
# the current row hits the cutoff, it resets the elements to the minimum value.
def generate_level_vectors(vector_length, cutoff):
    row = [-cutoff for i in range(0, vector_length)]
    ret = [row]

    # pops first iteration off the stack [-cutoff,-cutoff,-cutoff] which 
    # doesn't have the last value set and gives that to setMin function
    # that sets the correct last value for the permutation.  Note: this
    # only happens for the first permutation.
    first_permutation = ret.pop(0)
    setLastValue(first_permutation)
    print first_permutation
    while row[0:-1] != [cutoff for i in range(0, vector_length - 1)]:
        row = get_next_row(row, cutoff)
        ret.append(row)
    return ret

