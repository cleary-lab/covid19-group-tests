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
	args,_ = parser.parse_known_args()
	f = open(args.viral_load_matrix)
	ViralLoad = [[float(l) for l in line.strip().split()] for line in f]
	f.close()
	ViralLoad = 10**np.array(ViralLoad) - 1
	e0 = np.average([np.average(v[v>0] > args.LOD) for v in ViralLoad[:,args.start_time:args.end_time].T])
	Results = []
	FP = glob.glob(os.path.join(args.resultspath, 'Recall_combined.*.npy'))
	for fp in FP:
		cond = fp.split('.')[-2]
		n = int(cond.split('_')[0].split('-')[1])
		m = int(cond.split('_')[1].split('-')[1])
		sensitivity = np.load('%s/Recall_combined.%s.npy' % (args.resultspath,cond))
		efficiency = np.load('%s/Eff_avg.%s.npy' % (args.resultspath,cond))
		avg_sensitivity = np.average(sensitivity[args.start_time:args.end_time])
		avg_efficiency = np.average(efficiency[args.start_time:args.end_time])
		effective_m = n/avg_efficiency # stage (i) pools + stage (ii) validation tests
		Results.append((n,m,effective_m,avg_sensitivity))
	unique_nm = np.unique([r[0] for r in Results] + [r[1] for r in Results])
	X = np.zeros((len(unique_nm),len(unique_nm)))
	fig = plt.figure(num=None, figsize=(8, 8))
	for i,n_swabs in enumerate(unique_nm):
		for j,m_kits in enumerate(unique_nm):
			nm_min = min(n_swabs,m_kits)
			best = (e0*nm_min,nm_min,nm_min)
			for n,m,m_eff,sens in Results:
				n_tests = min(n_swabs/n, m_kits/m_eff)
				if n_tests > 0.9:
					effective_tests = n*n_tests*sens
					if effective_tests > best[0]:
						best = (effective_tests,n,m)
			X[j,i] = best[0]
			_=plt.text(i,j,'%d,%d' % (best[1],best[2]),horizontalalignment='center',verticalalignment='center',color="w",fontsize=8)
	_=plt.imshow(X,origin='lower',norm=mpl.colors.LogNorm())
	_=plt.xticks(np.arange(len(unique_nm)), labels=unique_nm)
	_=plt.yticks(np.arange(len(unique_nm)), labels=unique_nm)
	_=plt.xlabel('Number of swabs')
	_=plt.ylabel('Number of kits')
	_=plt.colorbar()
	_=plt.grid(b=None)
	_=plt.title('Effective individuals screened (n * sensitivity)')
	_=plt.tight_layout()
	_=plt.savefig('%s/summary.resource.png' % args.resultspath)
	plt.close()
		