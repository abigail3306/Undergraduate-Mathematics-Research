from math import floor
from math import ceil
from copy import copy
from generate_vectors import generate_level_vectors
from level_vectors import get_gap_vector
from level_vectors import one_line
from level_vectors import coxeter_length
import sys

################################
#  Utility functions
################################

# Given the base window w and an integer i,
# returns the image of i under the affine permutation w.
def value(w, i):
    ii = i % len(w)
    j = floor(i / len(w))
    if ii == 0:
        j = j-1
    # Note:  arrays in Python are indexed from zero, but we use the index 1 for the first element of our affine permutation.
    return w[ii-1] + j*len(w)

# Given a_list of values,
# returns the equivalent permutation in the same relative order
# containing the values 1, 2, ..., len(a_list).
def flatten(a_list):
    al_sorted = copy(a_list)
    al_sorted.sort()
    ret = []
    for i in range(0,len(a_list)):
        ret.append( al_sorted.index(a_list[i])+1 )
    return ret

# Given (the base window for) an affine permutation w,
# and an integer generator,
# returns the new affine permutation after acting by s_(generator).
# (for s_1, ..., s_{n-1}, we swap elements in positions i and i+1;
#  for s_0, we swap elements "across windows".)
def swap(w, generator):
    rlist = copy(w)
    size_of_window = len(rlist)
    intermediate_value = rlist[generator]
    if generator == 0:
        rlist[generator] = rlist[size_of_window - 1] - size_of_window
        rlist[size_of_window - 1] = intermediate_value + size_of_window
    else:
        rlist[generator] = rlist[generator - 1]
        rlist[generator - 1] = intermediate_value
    return rlist


################################
#  Pattern placement
################################

# This is a recursive function that tries to find an instance of the finite permutation pattern p in the affine permutation w.
# It places one entry from p during each pass.
# The placed parameter is the list of indices in w that represent the first positions in p.
# The a_vector_for_p records stopping conditions for each left-to-right maximum in p.  (It is just computed once and then passed without modification.)
def place_indecomposable(w, p, a_vector_for_p, placed):
    if len(placed) == len(p):
        return True

    prior_max = max(p[0:len(placed)])
    to_try = []
    
    # p[len(placed)] is the next entry of the pattern to try and place.
    # If it is a left-to-right maximum, we:
    #   find pair (a,b) with p[a] > p[b], a to the left of current pos, b to the right of current pos
    # Then once everything in a window of w is larger than the value we chose for p[a], we must stop (because it would be impossible to ever place p[b]).
    # Note that this pair (a,b) must exist for each left-to-right maximum since w is indecomposable.
    # TODO:  think about which choices of (a,b) are optimal?  Is there another stopping condition?
    if p[len(placed)] > prior_max:
        a = a_vector_for_p[len(placed)]
        if a == -1:
            print "ERROR:  a vector failure."
            return False

        prior_max_position = p.index(prior_max)
        prior_max_value = value(w, placed[ prior_max_position ])
        # si is the location in w to play the role of the next pattern entry.  
        # We initialize it to just after the last placed entry.
        si = placed[-1]+1
        # As explained above, we stop if the entire next window is larger than value(w, placed[a])
        while min([value(w, si+j) for j in range(0,len(w))]) < value(w, placed[a]):
            if value(w, si) > prior_max_value:
                to_try.append(si)
            si = si + 1

    # If p[len(placed)] is not a left-to-right maximum, then we can get bounds on where to stop looking in w directly from the pattern entries we've already placed.
    else:
        p_flat = flatten(p[0:len(placed)+1])
        upper_bound_p_position = p_flat.index( p_flat[len(placed)]+ 1)
        upper_bound_value = value(w, placed[upper_bound_p_position])
        if p_flat[len(placed)] > 1:
            lower_bound_p_position = p_flat.index( p_flat[len(placed)]- 1)
            lower_bound_value = value(w, placed[lower_bound_p_position])
        else:
            lower_bound_value = None

        ## We stop if the current value in w doesn't fit upper (and lower, if available) bounds from p.
        si = placed[-1]+1
        while min([value(w, si+j) for j in range(0,len(w))]) < upper_bound_value:
            if value(w, si) < upper_bound_value and (lower_bound_value == None or value(w, si) > lower_bound_value):
                to_try.append(si)
            si = si+1

    # In to_try, we've recorded all of the positions in w that are compatible with the next pattern entry.
    # We recursively try each of them.
    for j in to_try:
        extend_placed = copy(placed)
        extend_placed.append( j )
        if place_indecomposable(w, p, a_vector_for_p, extend_placed):
            return True
    return False

# This starts the recursive search for an instance of the finite pattern p inside the affine permutation w. 
# Since w is affine, we can assume without loss of generality that the first entry of p lies in the first window of w.
def contains_pattern(w, p, a_vector_for_p):
    for j in range(1,len(w)+1):
        if place_indecomposable(w, p, a_vector_for_p, [j]):
            return True
    return False


################################
#  Main entry point
################################

# TODO:  add an option to count only the sorted affine permutations (correspond to abaci).
# TODO:  can we count k-strand affine permutations?
# TODO:  can we parallelize this search?
# TODO:  conversely, can we make place_indecomposible somehow use the information earlier in the tree for speedup?
# TODO:  can we cache the full set of affine permutations and then run different patterns against them?

# TODO:  consolidate pre-processing to get bounds in both LR-max and non-LR-max cases.


# This computes the "direct sum" decomposition of the finite permutation p.
# p = a + b if all values in a lie below all the values in b.
# We can assume that p is indecomposible in our pattern placement algorithm because of the following:
#   Proposition:  Suppose p = a + b.  Then w contains p if and only if w contains a and b.
def sum_decompose(p):
    breakpos = [0]
    for i in range(0,len(p)):
        if set(p[0:i+1]) == set(range(1,i+2)):
            breakpos.append(i+1)
    return [ flatten(p[breakpos[i]:breakpos[i+1]]) for i in range(0,len(breakpos)-1) ]

# This is the main entry point for the code, called from "fastscope.py".
# It generates all affine permutations on N generators up to length k
# and keeps a count of those that avoid the finite permutation pattern p (if p is passed).
# Trace with p = 34251.  check indecomposible.  LR maxima at positions 1, 2, 4.  (optimal) a-values: n/a, 1, 3, respectively.
def affine_permutations_of_length(N, k, p=None):
    if p != None and len(sum_decompose(p)) > 1:
        print "Only indecomposable patterns are accepted."
        return

    # Precompute the a-vector for LR maxima
    # for each LR max, find pair (a,b) with p[a] > p[b], a to the left of current pos, b to the right of current pos
    # record a in the vector.
    prior_max = p[0]
    a_vector = [-1]
    for i in range(1,len(p)):
        # if the position is a left-to-right maximum
        if p[i] > prior_max:
            a = 0
            # TODO:  this is the first a value that works, not necessarily the minimal one.
            while a < i and len([ j for j in range(i+1 , len(p)) if p[j] < p[a] ]) == 0:
                a = a+1
            if a == i:
                print "ERROR:  p is decomposible: ", p, " at pos ", i
                return False
            a_vector.append(a)
            prior_max = p[i]
        # if the position is not a LR maximum, just record -1.
        else:
            a_vector.append(-1)

    if k == 0:
        return [[range(1,N+1)], 1]


####NEW IMPLEMENTATION#####
    # TODO: is there a way to determine minimal coxeter length (worst case?)
    # TODO: for any l(i) >= B, what is the lower bound for the coxeter length?
    
    # coxeter_list is the former p_avoiding_count.  I found it better to call
    # it this since it's no longer a single value but rather a list.
    coxeter_list = [0 for i in range(0, k+1)]

    # TODO: consider how to slim down list we're considering...
    # coxeter_sum = sum[0 for i in range(0,len(coxeter_list)]
    p_avoiding_count = 0
    running_total = 0
    cutoff = int(ceil(k/N)+1)
    affine_permutations = generate_level_vectors(N, cutoff)
   
    for x in affine_permutations:
	   running_total = running_total + 1
	   if not contains_pattern(one_line(x), p, a_vector):
	      print x
	      if coxeter_length(get_gap_vector(x)) < len(coxeter_list):
	         coxeter_list[coxeter_length(get_gap_vector(x))] = coxeter_list[coxeter_length(get_gap_vector(x))] + 1	  	    
	      p_avoiding_count = p_avoiding_count + 1
    
    print p_avoiding_count, " of length ", k, " avoiding ", p, ", ", running_total, "total."
    print coxeter_list    
