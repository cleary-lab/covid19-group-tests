import numpy as np
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
	parser.add_argument('--infection-time', help='Path to time of infection for each individual')
	parser.add_argument('--onset-time', help='Path to time of symptom onset for each individual')
	parser.add_argument('--resultspath', help='Path to save summary figure')
	parser.add_argument('--n-individuals', help='Number of individuals to test (n)',type=int)
	parser.add_argument('--m-pools', help='Number of pools (m)',type=int)
	parser.add_argument('--q-split', help='Number of pools (q) into which each sample is split',type=int)
	args,_ = parser.parse_known_args()
	f = open(args.viral_load_matrix)
	timepoints = f.readline()
	ViralLoad = [[float(l) for l in line.strip().split(',')] for line in f]
	f.close()
	ViralLoad = np.array(ViralLoad)
	x = [np.average(v[v > 0] > 2) for v in ViralLoad.T]
	f = open(args.infection_time)
	header = f.readline()
	InfectionTime = [float(line.strip()) for line in f]
	f.close()
	InfectionTime = np.array(InfectionTime)
	f = open(args.onset_time)
	header = f.readline()
	OnsetTime = [float(line.strip()) for line in f]
	f.close()
	OnsetTime = np.array(OnsetTime)
	cond = 'n-%d_m-%d_q-%d' % (args.n_individuals,args.m_pools,args.q_split)
	recall = np.load('%s/Recall_combined.%s.npy' % (args.resultspath,cond))
	efficiency = np.load('%s/Eff_avg.%s.npy' % (args.resultspath,cond))
	recall_n = np.load('%s/Recall_n.%s.npy' % (args.resultspath,cond),allow_pickle=True)[:-10] # exclude last 10d when infected rate > 1%
	Vl = np.load('%s/Vl.%s.npy' % (args.resultspath,cond),allow_pickle=True)[:-10] # exclude last 10d when infected rate > 1%
	Vi = np.load('%s/Vi.%s.npy' % (args.resultspath,cond),allow_pickle=True)[:-10] # exclude last 10d when infected rate > 1%
	Vo = np.load('%s/Vo.%s.npy' % (args.resultspath,cond),allow_pickle=True)[:-10] # exclude last 10d when infected rate > 1%
	v1,r1 = bin_results(np.hstack(recall_n), np.hstack(Vl))
	v2,r2 = bin_results(np.hstack(Vl), np.hstack(Vo),y0=100)
	v3,r3 = bin_results(np.hstack(recall_n), np.hstack(Vo))
	fig = plt.figure(num=None, figsize=(8, 10))
	ax1 = plt.subplot(3,2,1)
	_=plt.plot(recall,label='group test')
	_=plt.plot(x,label='individual test')
	_=plt.legend(loc='best',prop={'size': 8})
	_=plt.ylabel('Average sensitivity',fontsize=10)
	ax2 = plt.subplot(3,2,2)
	_=plt.plot(np.log10(v1),r1)
	_=plt.xlabel('log10( Viral load )',fontsize=10)
	_=plt.ylabel('Average sensitivity',fontsize=10)
	ax3 = plt.subplot(3,2,3)
	_=plt.plot(efficiency)
	_=plt.ylabel('Average efficiency',fontsize=10)
	ax4 = plt.subplot(3,2,4)
	_=plt.plot(v3,r3,label='group test')
	_=plt.plot(v2,1-r2,label='individual test')
	_=plt.legend(loc='best',prop={'size': 8})
	_=plt.xlabel('Time since symptom onset',fontsize=10)
	_=plt.ylabel('Average sensitivity',fontsize=10)
	ax5 = plt.subplot(3,2,5)
	_=plt.plot(efficiency*recall,label='group test')
	_=plt.plot(x,label='individual test')
	_=plt.legend(loc='best',prop={'size': 8})
	_=plt.xlabel('Time',fontsize=10)
	_=plt.ylabel('Total recall (efficiency*sensitivity)',fontsize=9)
	_=plt.tight_layout()
	plt.savefig('%s/summary.%s.pdf' % (args.resultspath,cond))
	plt.close()


