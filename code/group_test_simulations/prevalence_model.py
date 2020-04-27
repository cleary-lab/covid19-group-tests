import numpy as np
import argparse
import pandas as pd
from scipy.stats import poisson
import joblib
from estimate_prevalence import run_em
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns

def sample_pooled_ct_values(viral_load,batch_size,y_cept,Ct_thresh=40,LOD=1,fp_rate=0.01):
	idx = np.random.choice(viral_load.shape[0],batch_size,replace=False)
	vl = 10**viral_load[idx] - 1
	sampled_load = poisson.rvs(vl/batch_size)
	total_sampled_load = sampled_load.sum()
	if np.random.random() < fp_rate:
		# 1% of the time (or whatever fp_rate is) add the viral load from 1 positive sample selected at random with load sufficient to cause a false positive
		nz = np.where(viral_load > np.log10(LOD*batch_size))[0]
		i = np.random.choice(nz,1)
		total_sampled_load += poisson.rvs((10**viral_load[i]-1)/batch_size)
	ct_value = y_cept - np.log2(10)*np.log10(total_sampled_load+1e-8)
	if ct_value > Ct_thresh:
		ct_value = -1 # undetected
	sampled_freq = np.average(vl > LOD)
	return sampled_freq,ct_value

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--kde-path', help='Path to KDEs fit for each value of k, including prefix (ie, seir_kde)')
	parser.add_argument('--batch-size', help='Number of samples in each batch',type=int)
	parser.add_argument('--num-batches', help='Number of batches',type=int)
	parser.add_argument('--start-time', help='Start of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	parser.add_argument('--savepath', help='Path to save summary figure and results')
	parser.add_argument('--ct-intercept-mean', help='Mean of Ct y-intercept distribution',type=float,default=38.5)
	parser.add_argument('--ct-intercept-std', help='Stdev of Ct y-intercept distribution',type=float,default=1)
	parser.add_argument('--fp-rate', help='Rate of false positive PCRs',type=float,default=0.01)
	parser.add_argument('--fp-rate-estimate', help='Estimated false positive PCR rate',type=float,default=0.002)
	parser.add_argument('--p0', help='Initial frequency estimate',type=float,default=0.005)
	parser.add_argument('--num-trials', help='Number of random sample trials at each time point',type=int,default=100)
	parser.add_argument('--make-plot',dest='make_plot', action='store_true')
	parser.set_defaults(make_plot=False)
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	N = args.batch_size*args.num_batches
	ViralLoad_test = np.load(args.viral_load_matrix)
	Kernels = [joblib.load('%s.k-%d.pkl' % (args.kde_path,k)) for k in range(1,args.batch_size+1)]
	Results = []
	true_frequencies = []
	for ti in range(args.start_time,args.end_time):
		p = []
		for _ in range(args.num_trials):
			p_obs = []
			Y = []
			y_cept = np.random.randn()*args.ct_intercept_std + args.ct_intercept_mean
			for b in range(args.num_batches):
				freq,y = sample_pooled_ct_values(ViralLoad_test[:,ti],args.batch_size,y_cept,fp_rate=args.fp_rate)
				p_obs.append(freq)
				Y.append(y)
			p_obs = np.average(p_obs)
			if max(Y) > 0:
				p_est = run_em(Y,args.num_batches,args.batch_size,args.p0,Kernels,args.fp_rate_estimate)
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

