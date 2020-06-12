import numpy as np
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--savepath', help='Path to save grid params')
	parser.add_argument('--n-individuals', help='Number of individuals to test (n)',type=int)
	args,_ = parser.parse_known_args()
	
	# find all factors
	n = args.n_individuals
	factors = np.flatnonzero(n % np.arange(1,n) == 0)+1
	
	# Compute unique pairs
	pairs = set(tuple(sorted([f,n//f])) for f in factors)
	
	# Append to file
	with open(args.savepath, 'a') as f:
		for p in pairs:
			f.write('%s\t%s\n' % p)

# To generate, run the following:
# python -u gengridsweep.py --n-individuals 96 --savepath grid_params.txt
# python -u gengridsweep.py --n-individuals 192 --savepath grid_params.txt