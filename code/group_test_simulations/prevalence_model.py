import numpy as np
import argparse
import pandas as pd
from scipy.special import comb
from scipy.stats import poisson
from sklearn.neighbors import KernelDensity
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns

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

def total_density(m,p,conditional_k):
	x = 0
	for k in range(m+1):
		x += comb(m,k)*p**k*(1-p)**(m-k)*conditional_k[k]
	return x

def cond_prob(p,m,k,td,conditional_k):
	cp = comb(m,k)*p**k*(1-p)**(m-k)*conditional_k[k]/td
	return cp

def run_em(y,b,n,p0,Kernels,max_itr=100,min_itr=10,tol=1e-3,n_vals=50):
	# calculate all the conditional densities for the observed values (y) up front
	CK = []
	ck0 = [1]+[np.exp(Kernels[k].score_samples([[-1]])[0]) for k in range(n)] # conditional density of 0.1
	for i in range(b):
		if y[i] > 0:
			conditional_k = [0]+[np.exp(Kernels[k].score_samples([[np.log10(y[i]*n)]])[0]) for k in range(n)] # conditional density of y*m, given k
		else:
			conditional_k = ck0
		CK.append(conditional_k)
	for itr in range(max_itr):
		p = 0
		for i in range(b):
			conditional_k = CK[i]
			td = total_density(n,p0,conditional_k)
			for k in range(1,n+1):
				p += k*cond_prob(p0,n,k,td,conditional_k)
		p = p/(n*b)
		if (itr > min_itr) and (abs(p-p0)/p0 < tol):
			break
		p0 = p
	return p

def sample_pooled_viral_load(viral_load,batch_size,LOD=1):
	idx = np.random.choice(viral_load.shape[0],batch_size,replace=False)
	vl = 10**viral_load[idx] - 1
	sampled_load = poisson.rvs(vl/batch_size)
	total_sampled_load = sampled_load.sum()
	if total_sampled_load < LOD:
		total_sampled_load = 0
	sampled_freq = np.average(vl > LOD)
	return sampled_freq,total_sampled_load

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--batch-size', help='Number of samples in each batch',type=int)
	parser.add_argument('--num-batches', help='Number of batches',type=int)
	parser.add_argument('--start-time', help='Start of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	parser.add_argument('--savepath', help='Path to save summary figure and results')
	parser.add_argument('--p0', help='Initial frequency estimate',type=float,default=0.005)
	parser.add_argument('--num-trials', help='Number of random sample trials at each time point',type=int,default=100)
	parser.add_argument('--make-plot',dest='make_plot', action='store_true')
	parser.set_defaults(make_plot=False)
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	N = args.batch_size*args.num_batches
	f = open(args.viral_load_matrix)
	timepoints = f.readline()
	ViralLoad = [[float(l) for l in line.strip().split(',')] for line in f]
	f.close()
	ViralLoad = np.array(ViralLoad)
	print('matrix size: %d x %d' % (ViralLoad.shape[0],ViralLoad.shape[1]))
	# use 20% of the data for learning empirical distributions
	train_idx = np.random.choice(ViralLoad.shape[0],int(ViralLoad.shape[0]/5),replace=False)
	test_idx = np.setdiff1d(np.arange(ViralLoad.shape[0]),train_idx)
	ViralLoad_train = ViralLoad[train_idx]
	ViralLoad_test = ViralLoad[test_idx]
	Kernels = empirical_convolutions(ViralLoad_train,args.batch_size)
	Results = []
	true_frequencies = []
	for ti in range(args.start_time,args.end_time):
		p = []
		for _ in range(args.num_trials):
			p_obs = []
			Y = []
			for b in range(args.num_batches):
				freq,y = sample_pooled_viral_load(ViralLoad_test[:,ti],args.batch_size)
				p_obs.append(freq)
				Y.append(y)
			p_obs = np.average(p_obs)
			if max(Y) > 0:
				p_est = run_em(Y,args.num_batches,args.batch_size,args.p0,Kernels)
			else:
				p_est = 0
			p.append([p_obs,p_est])
		p = np.array(p)
		Results.append(p)
		true_freq = np.average(ViralLoad_test[:,ti] > 0)
		true_frequencies.append(true_freq)
		print('Time point: %d; Population frequency: %.5f; Sampled frequency: %.5f; Estimated frequency: %.5f' % (ti, true_freq, np.average(p[:,0]), np.average(p[:,1])))
	Results = np.array(Results) + 1e-5 # we'll say the baseline estimated frequency is 1/100,000
	np.save('%s/results.b-%d_n-%d_t0-%d_t1-%d.npy' % (args.savepath,args.num_batches,args.batch_size,args.start_time,args.end_time),Results)
	np.save('%s/true_frequencies.b-%d_n-%d_t0-%d_t1-%d.npy' % (args.savepath,args.num_batches,args.batch_size,args.start_time,args.end_time),true_frequencies)
	if args.make_plot:
		mode = np.empty(Results.shape,dtype='<U42')
		mode[:,:,0] = 'sampled population; N=%d,b=%d,n=%d' % (N,args.num_batches,args.batch_size)
		mode[:,:,1] = 'estimated from b=%d measurements' % args.num_batches
		F = (np.ones(Results.shape).T*np.array(true_frequencies)).T
		df = {'true frequency of positive samples': F.flatten(), 'estimated frequency of positive samples': Results.flatten(), 'mode': mode.flatten()}
		df = pd.DataFrame(df)
		ax = sns.lineplot(data=df, x='true frequency of positive samples', y='estimated frequency of positive samples', hue='mode', ci='sd')
		ax.set_xscale('log')
		ax.set_yscale('log')
		handles, labels = ax.get_legend_handles_labels()
		_= ax.legend(handles=handles[1:], labels=labels[1:]) # remove legend title
		_=plt.title('True vs Estimated prevalence, analyzed from days %d-%d' % (args.start_time,args.end_time))
		_=plt.tight_layout()
		plt.savefig('%s/estimated_prevalence.b-%d_n-%d_t0-%d_t1-%d.pdf' % (args.savepath,args.num_batches,args.batch_size,args.start_time,args.end_time))
		plt.close()

