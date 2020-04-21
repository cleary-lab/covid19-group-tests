import numpy as np
import argparse
from sklearn.neighbors import KernelDensity
from sklearn.externals import joblib

def empirical_convolutions(ViralLoad,K,b=0.1,samples=10000,thresh=0):
	Kernels = []
	nz = np.where(ViralLoad.flatten() > thresh)[0]
	x = ViralLoad.flatten()[nz]
	x = 10**x
	for k in range(1,K+1):
		z = [np.random.choice(x,k).sum() for _ in range(samples)]
		z = np.log10(z)
		kernel = KernelDensity(bandwidth=b).fit(z[:,np.newaxis])
		Kernels.append(kernel)
	return Kernels

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--batch-size', help='Number of samples in each batch',type=int)
	parser.add_argument('--savepath', help='Path to save summary figure and results')
	parser.add_argument('--bandwidth', help='KDE bandwidth',type=float,default=0.1)
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	f = open(args.viral_load_matrix)
	timepoints = f.readline()
	ViralLoad = [[float(l) for l in line.strip().split(',')] for line in f]
	f.close()
	ViralLoad = np.array(ViralLoad)
	Kernels = empirical_convolutions(ViralLoad,args.batch_size,b=args.bandwidth)
	for k in range(1,args.batch_size+1):
		_=joblib.dump(Kernels[k-1], '%s/seir_kde.k-%d.pkl' % (args.savepath,k))
