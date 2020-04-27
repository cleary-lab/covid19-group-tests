import numpy as np
import argparse
from sklearn.neighbors import KernelDensity
from sklearn.externals import joblib

def get_random_intercepts(mean_range,std_range,n):
	# normal distribution with uniform prior on mean and std
	means = np.random.random(n)*(mean_range[1]-mean_range[0]) + mean_range[0]
	stds = np.random.random(n)*(std_range[1]-std_range[0]) + std_range[0]
	y_cepts = np.random.randn(n)*stds + means
	return y_cepts

def empirical_convolutions(ViralLoad,K,b=0.1,samples=10000,thresh=0,mean_range=(37,39),std_range=(0.1,1)):
	Kernels = []
	nz = np.where(ViralLoad.flatten() > thresh)[0]
	x = ViralLoad.flatten()[nz]
	x = 10**x
	for k in range(1,K+1):
		z = [np.random.choice(x,k).sum() for _ in range(samples)]
		z = np.log10(z)
		# convert log10( viral load ) to Ct value
		y_cepts = get_random_intercepts(mean_range,std_range,samples)
		ct_values = y_cepts - np.log2(10)*z
		ct_values[ct_values < 0] = 0
		kernel = KernelDensity(bandwidth=b).fit(ct_values[:,np.newaxis])
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
	# use 20% of the data for learning empirical distributions
	train_idx = np.random.choice(ViralLoad.shape[0],int(ViralLoad.shape[0]/5),replace=False)
	test_idx = np.setdiff1d(np.arange(ViralLoad.shape[0]),train_idx)
	ViralLoad_train = ViralLoad[train_idx]
	ViralLoad_test = ViralLoad[test_idx]
	Kernels = empirical_convolutions(ViralLoad_train,args.batch_size,b=args.bandwidth)
	for k in range(1,args.batch_size+1):
		_=joblib.dump(Kernels[k-1], '%s/seir_kde.k-%d.pkl' % (args.savepath,k))
	np.save('%s.train.npy' % (args.viral_load_matrix[:args.viral_load_matrix.rfind('.')]), ViralLoad_train)
	np.save('%s.test.npy' % (args.viral_load_matrix[:args.viral_load_matrix.rfind('.')]), ViralLoad_test)
