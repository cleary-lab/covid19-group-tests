import numpy as np
import glob,os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns
import argparse

def bin_results(R,P0,n=50,y0=None):
	P = np.copy(P0)
	lower,upper = np.percentile(P,1),np.percentile(P,99)
	P[P<lower] = lower
	P[P>upper] = upper
	xs = np.argsort(P)
	bins = [P[xs[i]] for i in range(0,len(P),int(len(P)/n))]
	indices = np.digitize(P,bins)
	x = np.array([np.average(P[indices == i]) for i in set(indices)])
	if y0 is not None:
		y = np.array([np.average(R[indices == i] < y0) for i in set(indices)])
	else:
		y = np.array([np.average(R[indices == i]) for i in set(indices)])
	return x,y

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--resultspath', help='Path to save summary figure')
	parser.add_argument('--start-time', help='Beginning of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	parser.add_argument('--LOD', help='Limit of detection',type=float,default=100)
	parser.add_argument('--min-resources', help='Minimum number of resources to consider',type=int,default=12)
	args,_ = parser.parse_known_args()
	f = open(args.viral_load_matrix)
	ViralLoad = [[float(l) for l in line.strip().split()] for line in f]
	f.close()
	ViralLoad = 10**np.array(ViralLoad) - 1
	e0 = np.average([np.average(v[v>0] > args.LOD) for v in ViralLoad[:,args.start_time:args.end_time].T if v.sum()!=0])
	Results = []
	unique_nm = []
	FP = glob.glob(os.path.join(args.resultspath, 'Recall_combined.*.npy'))
	for fp in FP:
		cond = fp.split('.')[-2]
		n = int(cond.split('_')[0].split('-')[1])
		m = int(cond.split('_')[1].split('-')[1])
		q = int(cond.split('_')[2].split('-')[1])
		q = min(m,q) # filenames sometime screwy
		sensitivity = np.load('%s/Recall_combined.%s.npy' % (args.resultspath,cond))
		efficiency = np.load('%s/Eff_avg.%s.npy' % (args.resultspath,cond))
		avg_sensitivity = sensitivity[args.start_time:args.end_time]
		avg_efficiency = efficiency[args.start_time:args.end_time]
		#effective_m = n/avg_efficiency # stage (i) pools + stage (ii) validation tests
		Results.append((n,m,q,avg_sensitivity,avg_efficiency))
		#print(Results[-1],m/Results[-1][3]*Results[-1][4])
		if (n >= args.min_resources) and (m >= args.min_resources):
			unique_nm += [n,m]
	unique_nm = np.unique(unique_nm)
	X = np.zeros((len(unique_nm),len(unique_nm)))
	fig = plt.figure(num=None, figsize=(8, 8))
	for i,n_swabs in enumerate(unique_nm):
		for j,m_kits in enumerate(unique_nm):
			nm_min = min(n_swabs,m_kits)
			best = (e0*nm_min,nm_min,nm_min,1,1)
			for n,m,q,sens,eff in Results:
				m_eff = n/eff # stage (i) pools + stage (ii) validation tests
				n_tests = np.minimum(n_swabs/n, m_kits/m_eff)
				if np.average(n_tests) > 0.9:
					effective_tests = np.average(n*n_tests*sens)
					if effective_tests > best[0]:
						best = (effective_tests,n,m,q,np.average(n_tests))
			X[j,i] = best[0]/(e0*nm_min)
			if X[j,i] < 10:
				# white font
				_=plt.text(i,j,'%d\n%d\n%d\n%.1f' % (best[1],best[2],best[3],best[4]), horizontalalignment='center',verticalalignment='center',color="w",fontsize=7)
			else:
				# black font
				_=plt.text(i,j,'%d\n%d\n%d\n%.1f' % (best[1],best[2],best[3],best[4]), horizontalalignment='center',verticalalignment='center',color="black",fontsize=7)
	print(X.min(),X.max())
	_=plt.imshow(X,origin='lower',cmap='cubehelix')
	_=plt.xticks(np.arange(len(unique_nm)), labels=unique_nm)
	_=plt.yticks(np.arange(len(unique_nm)), labels=unique_nm)
	_=plt.xlabel('Daily samples collected')
	_=plt.ylabel('Daily reaction (RNA extraction, RT+qPCR) budget')
	_=plt.colorbar()
	_=plt.grid(b=None)
	_=plt.title('Effectiveness of optimal design from days %d-%d' % (args.start_time,args.end_time))
	_=plt.tight_layout()
	_=plt.savefig('%s/summary.resource.pdf' % args.resultspath)
	plt.close()
		