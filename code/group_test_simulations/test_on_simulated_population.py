import numpy as np
import spams
from scipy.stats import poisson
import glob,os
import argparse

def random_binary_balanced(m,n,q):
	A = np.zeros((m,n))
	for i in range(n):
		idx = np.random.choice(m,q,replace=False)
		A[idx,i] = 1
	return A

def decode(y,A,e):
	x = A.sum(0) - A.T.dot(y)
	return (x < e)

THREADS=40
def sparse_decode(Y,D,lda):
	Ynorm = np.linalg.norm(Y)**2/Y.shape[1]
	W = spams.lasso(np.asfortranarray(Y),np.asfortranarray(D),lambda1=lda*Ynorm,mode=1,numThreads=THREADS,pos=True)
	W = np.asarray(W.todense())
	return W

def determine_sampling_volume(A,max_volume):
	max_samples_per_pool = A.sum(1).max()
	return np.round(max_volume/max_samples_per_pool,1)

def sample_pooled_viral_load(A_i,viral_load,vol_total,vol_sampled):
	sampled_load = poisson.rvs(viral_load*vol_sampled/vol_total)
	total_sampled_load = A_i.dot(sampled_load)
	return total_sampled_load

def get_n_random_samples(A,ViralLoad,vol_total,vol_sampled,LOD):
	# should consider not sampling individuals who have had past infections
	idx_n = np.random.choice(ViralLoad.shape[0],A.shape[1],replace=False)
	viral_load = ViralLoad[idx_n]
	x_true = (viral_load > 0).astype(np.float) # should this be > LOD?
	y = np.zeros(A.shape[0])
	for i in range(len(y)):
		total_sampled_load = sample_pooled_viral_load(A[i],viral_load,vol_total,vol_sampled)
		if total_sampled_load > LOD:
			y[i] = 1
	return x_true,y,idx_n

def run_test(A,ViralLoad,ViralTime,vol_total,vol_sampled,LOD,e,itr=2500):
	true_pos = 0
	total_pos = 0
	T = []
	Vl = []
	Vt = []
	R = []
	for __ in range(itr):
		x,y,idx_n = get_n_random_samples(A,ViralLoad,vol_total,vol_sampled,LOD)
		# we take the union of a naive decoder (decode) and a lasso (sparse_decode)
		# the naive decoder is doing the bulk of the work here, but occasionally the lasso picks up a few extra hits
		x_hat1 = decode(y,A,e)
		x_hat2 = sparse_decode(y[:,np.newaxis],A,0.01)
		x_hat2 = (x_hat2[:,0] > 0.01)
		x_hat3 = ((x_hat1 + x_hat2) > 0)
		true_pos += x.dot(x_hat3)
		total_pos += x.sum()
		t = x_hat3.sum()
		T.append(t)
		if x.sum() > 0:
			ki = np.where(x)[0]
			R.append(x[ki]*x_hat3[ki])
			Vl.append(ViralLoad[idx_n][ki])
			Vt.append(ViralTime[idx_n][ki])
	T = np.array(T)
	if len(R):
		R = np.hstack(R)
		Vl = np.hstack(Vl)
		Vt = np.hstack(Vt)
	efficiency_avg = A.shape[1]/(A.shape[0] + np.average(T))
	if total_pos > 0:
		recall = true_pos/total_pos
	else:
		recall = 1
	return efficiency_avg, recall, R, Vl, Vt

def relative_infection_time(ViralLoad):
	# only writing the viral time after infection for the days with non-zero load here (all other days set to -1)
	time_of_infection = {}
	ViralTime = np.ones(ViralLoad.shape)*-1
	for i,t in zip(*np.where(ViralLoad)):
		if i not in time_of_infection:
			time_of_infection[i] = t
		ViralTime[i,t] = t-time_of_infection[i]
	return ViralTime

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--expected-pos-prob', help='1-expected faulty test rate (for setting error correction threshold)',type=float)
	parser.add_argument('--savepath', help='Path to save summary figure')
	parser.add_argument('--n-individuals', help='Number of individuals to test (n)',type=int)
	parser.add_argument('--m-pools', help='Number of pools (m)',type=int)
	parser.add_argument('--max-dilution', help='Maximum number of samples per composition',type=int,default=100)
	parser.add_argument('--max-sample-split', help='Maximum number of compositions a single sample can be split into',type=int,default=8)
	parser.add_argument('--max-pool-volume', help='Maximum volume of pooled samples',type=int,default=140)
	parser.add_argument('--volume-per-swab', help='Volume used per swab in individual test',type=int,default=100)
	parser.add_argument('--LOD', help='Limit of detection',type=float,default=100)
	args,_ = parser.parse_known_args()
	f = open(args.viral_load_matrix)
	ViralLoad = [[float(l) for l in line.strip().split()] for line in f]
	f.close()
	ViralLoad = 10**np.array(ViralLoad) - 1
	ViralTime = relative_infection_time(ViralLoad)
	E_avg = np.zeros(ViralLoad.shape[1])
	Recall_combined = np.zeros(ViralLoad.shape[1])
	Recall_n, Vl, Vt = [],[],[]
	q = int(min(args.max_dilution/args.n_individuals*args.m_pools,args.max_sample_split,args.m_pools/2,np.round(np.log(args.n_individuals))))
	# pooling composition matrix
	A = random_binary_balanced(args.m_pools,args.n_individuals,q)
	# number of faulty tests to tolerate
	e = np.round(2*(1-args.expected_pos_prob)*q) + 1e-7
	vol_sampled = determine_sampling_volume(A,args.max_pool_volume)
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	print('q\t%d' % q)
	print('e\t%d' % e)
	print('volume sampled\t%.1f' % vol_sampled)
	for t in range(ViralLoad.shape[1]):
		if ViralLoad[:,t].sum() > 0:
			ea,rec,r,vl,vt = run_test(A,ViralLoad[:,t],ViralTime[:,t],args.volume_per_swab,vol_sampled,args.LOD,e)
			Recall_n.append(r)
			Vl.append(vl)
			Vt.append(vt)
		else:
			# save time by not running simulation if there are no cases
			# will want to change this if we add technical (eg pcr) false positives
			ea,rec = args.n_individuals/args.m_pools, 1
		E_avg[t] = ea
		Recall_combined[t] = rec
		print('Time: %d; Case frequency: %.5f; Efficiency: %.2f; Recall: %.4f; Fraction >LOD: %.4f' % (t,np.average(ViralLoad[:,t] > 0),ea,rec, np.average(ViralLoad[ViralLoad[:,t] > 0,t] > args.LOD)))
	np.save('%s/Eff_avg.n-%d_m-%d.npy' % (args.savepath,args.n_individuals,args.m_pools),E_avg)
	np.save('%s/Recall_combined.n-%d_m-%d.npy' % (args.savepath,args.n_individuals,args.m_pools),Recall_combined)
	np.save('%s/Recall_n.n-%d_m-%d.npy' % (args.savepath,args.n_individuals,args.m_pools),Recall_n)
	np.save('%s/Vl.n-%d_m-%d.npy' % (args.savepath,args.n_individuals,args.m_pools),Vl)
	np.save('%s/Vt.n-%d_m-%d.npy' % (args.savepath,args.n_individuals,args.m_pools),Vt)

