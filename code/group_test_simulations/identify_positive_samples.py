import numpy as np
import argparse

def decode(y,A,e):
	x = A.sum(0) - A.T.dot(y)
	return (x < e)

def identify_putative_cases(y,A,e=0.1):
	x_hat1 = decode(y,A,e)
	return np.where(x_hat1)[0]

def parse_composition_file(fp):
	f = open(fp)
	A_idx = [[int(x)-1 for x in line.strip().split(',')] for line in f]
	f.close()
	num_samples = len(A_idx)
	num_pools = max(max(idx) for idx in A_idx) + 1
	A = np.zeros((num_pools,num_samples))
	for i,idx in enumerate(A_idx):
		A[idx,i] = 1
	return A

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--observations', help='Path to file with positive (1) or negative (0) result of each pooled test')
	parser.add_argument('--pools', help='Path to pool compositions. One line per sample with comma-separated list of pools')
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	f = open(args.observations)
	Y = np.array([float(line.strip()) for line in f])
	f.close()
	A = parse_composition_file(args.pools)
	putatives = identify_putative_cases(Y,A)
	print('\nFound %d putative positives:' % len(putatives))
	for i in putatives:
		print('sample %d' % (i+1))
	
	
	