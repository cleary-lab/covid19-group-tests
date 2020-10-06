import numpy as np
from scipy.sparse import save_npz, csr_matrix
import argparse

def parse_sparse_viral_loads(fp,N):
	rows = []
	cols = []
	vals = []
	f = open(fp)
	header = f.readline()
	for line in f:
		if len(line) < 100: # some corrupted lines with tons of characters
			i,t,vl = line.strip().split(',')
			rows.append(int(i)-2) # indices are offset by 1, plus one indexed
			cols.append(int(t))
			vals.append(float(vl))
	f.close()
	vals = 10**np.array(vals) - 1
	T = max(cols)+1
	timepoints = np.arange(T)
	X = csr_matrix((vals,(rows,cols)), shape=(N,T))
	peak_time = X.argmax(axis=1)[:,0]
	return X.tocsc(),timepoints,peak_time

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--infection-time', help='Path to time of infection for each individual')
	args,_ = parser.parse_known_args()
	f = open(args.infection_time)
	header = f.readline()
	InfectionTime = []
	for line in f:
		if '.' in line: # non-empty lines
			InfectionTime.append(float(line.strip().split(',')[-1]))
		else:
			InfectionTime.append(-1)
	f.close()
	print('parsed infection times for %d individuals' % len(InfectionTime))
	ViralLoad, timepoints, PeakTime = parse_sparse_viral_loads(args.viral_load_matrix, len(InfectionTime))
	print('parsed viral load matrix with shape (%d, %d) and %d entries' % (ViralLoad.shape[0], ViralLoad.shape[1], len(ViralLoad.data)))
	outpath = args.viral_load_matrix[:args.viral_load_matrix.rfind('.')]
	print('Saving to %s' % outpath)
	save_npz('%s.viral_loads.npz' % outpath, ViralLoad)
	np.save('%s.timepoints.npy' % outpath, timepoints)
	np.save('%s.peak_times.npy' % outpath, PeakTime)
