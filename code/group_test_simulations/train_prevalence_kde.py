import numpy as np
import argparse
from sklearn.neighbors import KernelDensity
from sklearn.externals import joblib
from scipy.sparse import load_npz, save_npz

def get_random_intercepts(mean_range,std_range,n):
	# normal distribution with uniform prior on mean and std
	means = np.random.random(n)*(mean_range[1]-mean_range[0]) + mean_range[0]
	stds = np.random.random(n)*(std_range[1]-std_range[0]) + std_range[0]
	y_cepts = np.random.randn(n)*stds + means
	return y_cepts

def softplus(x):
	return np.log(1 + np.exp(x))

def empirical_convolutions(ViralLoad,K,b=0.1,samples=10000,thresh=1,mean_range=(37,39),std_range=(0.1,1)):
	Kernels = []
	x = ViralLoad.data[ViralLoad.data > thresh]
	x[x > 1e16] = 1e16
	for k in range(1,K+1):
		z = [np.random.choice(x,k).sum() for _ in range(samples)]
		z = np.log10(z)
		# convert log10( viral load ) to Ct value
		y_cepts = get_random_intercepts(mean_range,std_range,samples)
		ct_values = y_cepts - np.log2(10)*z
		# soft plus
		ct_values = softplus(ct_values)
		kernel = KernelDensity(bandwidth=b).fit(ct_values[:,np.newaxis])
		Kernels.append(kernel)
	return Kernels

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--batch-size', help='Number of samples in each batch',type=int)
	parser.add_argument('--savepath', help='Path to save summary figure and results')
	parser.add_argument('--growth-end', help='End of growth period',type=int,default=130)
	parser.add_argument('--decline-start', help='Start of decline period',type=int,default=155)
	parser.add_argument('--bandwidth', help='KDE bandwidth',type=float,default=0.1)
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	ViralLoad = load_npz(args.viral_load_matrix)
	# use 20% of the data for learning empirical distributions
	train_idx = np.random.choice(ViralLoad.shape[0],int(ViralLoad.shape[0]/5),replace=False)
	test_idx = np.setdiff1d(np.arange(ViralLoad.shape[0]),train_idx)
	ViralLoad_train = ViralLoad[train_idx]
	ViralLoad_test = ViralLoad[test_idx]
	print('estimating growth convolutions')
	Kernels = empirical_convolutions(ViralLoad_train[:,:args.growth_end],args.batch_size,b=args.bandwidth)
	print('saving')
	for k in range(1,args.batch_size+1):
		_=joblib.dump(Kernels[k-1], '%s_growth/seir_kde.k-%d.pkl' % (args.savepath,k))
	print('estimating decline convolutions')
	Kernels = empirical_convolutions(ViralLoad_train[:,args.decline_start:],args.batch_size,b=args.bandwidth)
	print('saving')
	for k in range(1,args.batch_size+1):
		_=joblib.dump(Kernels[k-1], '%s_decline/seir_kde.k-%d.pkl' % (args.savepath,k))
	save_npz('%s.train.npz' % (args.viral_load_matrix[:args.viral_load_matrix.rfind('.')]), ViralLoad_train)
	save_npz('%s.test.npz' % (args.viral_load_matrix[:args.viral_load_matrix.rfind('.')]), ViralLoad_test)
