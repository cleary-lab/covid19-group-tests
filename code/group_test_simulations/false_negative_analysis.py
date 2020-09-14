import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--resultspath', help='Path to save summary figure')
	parser.add_argument('--n-individuals', help='Number of individuals to test (n)',type=int)
	parser.add_argument('--m-pools', help='Number of pools (m)',type=int)
	parser.add_argument('--q-split', help='Number of pools (q) into which each sample is split',type=int)
	args,_ = parser.parse_known_args()
	cond = 'n-%d_m-%d_q-%d' % (args.n_individuals,args.m_pools,args.q_split)
	recall = np.load('%s/Recall_combined.%s.npy' % (args.resultspath,cond))
	efficiency = np.load('%s/Eff_avg.%s.npy' % (args.resultspath,cond))
	R = np.load('%s/Recall_n.%s.npy' % (args.resultspath,cond),allow_pickle=True)
	Vl = np.load('%s/Vl.%s.npy' % (args.resultspath,cond),allow_pickle=True)
	Vi = np.load('%s/Vi.%s.npy' % (args.resultspath,cond),allow_pickle=True)
	Vo = np.load('%s/Vo.%s.npy' % (args.resultspath,cond),allow_pickle=True)
	### peak times were originally saved without time offset...remove this in the future
	print('WARNING!!')
	print('peak times were originally saved without time offset')
	print('offsetting by range(10,121)')
	print('remove this in the future')
	Vo = np.array([t-v for t,v in zip(range(10,121),Vo)])
	###
	idx_false_neg = [np.where((r == 0)*(v < 1e12))[0] for r,v in zip(R,Vl)]
	peak_times_neg = np.hstack([v[idx] for v,idx in zip(Vo,idx_false_neg)])
	viral_loads_neg = np.log10(np.hstack([v[idx] for v,idx in zip(Vl,idx_false_neg)]))
	idx_true_pos = [np.where((r == 1)*(v < 1e12))[0] for r,v in zip(R,Vl)]
	peak_times_pos = np.hstack([v[idx] for v,idx in zip(Vo,idx_true_pos)])
	viral_loads_pos = np.log10(np.hstack([v[idx] for v,idx in zip(Vl,idx_true_pos)]))
	print('Proportion of false negatives before peak viral load: %f' % np.average(peak_times_neg < 0))
	print('Proportion of false negatives during first week following peak viral load: %f' % np.average((peak_times_neg >= 0)*(peak_times_neg < 7)))
	print('Proportion of false negatives 7 or more days after peak viral load: %f' % np.average(peak_times_neg >= 7))
	g = sns.JointGrid(peak_times_neg, viral_loads_neg, xlim=(-4,25), ylim=(-2,10))
	ax = sns.kdeplot(peak_times_pos,viral_loads_pos, cmap="Blues", ax=g.ax_joint, bw=(1,0.25), levels=15,shade=True,shade_lowest=False)
	ax = sns.kdeplot(peak_times_neg,viral_loads_neg, cmap="Reds", ax=g.ax_joint, bw=(1,0.25), levels=15,shade=True,shade_lowest=False)
	ax.set(xlabel='Time since peak viral load (days)', ylabel='log10( viral copies )')
	kwargs = {'cumulative': False}
	ax = sns.distplot(peak_times_neg, bins=30, kde_kws=kwargs, color="r", ax=g.ax_marg_x)
	ax = sns.distplot(peak_times_pos, bins=30, kde_kws=kwargs, color="b", ax=g.ax_marg_x)
	ax = sns.distplot(viral_loads_neg, kde_kws=kwargs, color="r", ax=g.ax_marg_y, vertical=True)
	ax = sns.distplot(viral_loads_pos, kde_kws=kwargs, color="b", ax=g.ax_marg_y, vertical=True)
	_=plt.tight_layout()
	plt.savefig('%s/summary.false_negative_density.%s.pdf' % (args.resultspath,cond))
	plt.close()