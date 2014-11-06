from scope import affine_permutations_of_length
import sys

if len(sys.argv) == 3:
    affine_permutations_of_length(int(sys.argv[1]), int(sys.argv[2]))
elif len(sys.argv) == 4:
    affine_permutations_of_length(int(sys.argv[1]), int(sys.argv[2]), list(eval(sys.argv[3])))
else:
    print "Usage:  python scope.py <number of generators> <length> <pattern to avoid (optional)>"
                       
