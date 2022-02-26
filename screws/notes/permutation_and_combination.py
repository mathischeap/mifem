# A Python program to print all
# permutations of given length
from itertools import permutations

# Get all permutations of length 2
# and length 2
perm = permutations([1, 2, 3], 2)

# Print the obtained permutations
for i in list(perm):
    print(i)




print('\n\n')


# A Python program to print all
# combinations of a given length
from itertools import combinations

# Get all combinations of [1, 2, 3]
# and length 2
comb = combinations([1, 2, 3], 2)

# Print the obtained combinations
for i in list(comb):
    print(i)