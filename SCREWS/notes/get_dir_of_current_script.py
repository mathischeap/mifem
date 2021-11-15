






# way 1
import os
current_dir = os.path.dirname(__file__)
print(current_dir.__class__.__name__, current_dir)


# way 2
import sys
current_dir = os.path.dirname(sys.argv[0])
print(current_dir.__class__.__name__, current_dir)