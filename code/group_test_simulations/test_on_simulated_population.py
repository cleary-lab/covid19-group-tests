import numpy as np
from scipy.stats import poisson
import glob,os
import argparse

# based on colbourn2006hoc VII.5.46
def genpairs(n):
	# Create starter design
	# starter = [(0,2*n-1)] + [(i,-i % (2*n-1)) for i in range(1,n)]
	Astart = np.zeros((2*n,n))
	Astart[[0,2*n-1],0] = 1
	for i in range(1,n):
		Astart[[i,-i % (2*n-1)],i] = 1
	
	# Add cycled versions
	A = np.copy(Astart)
	for r in range(1,2*n-1):
		Astart = np.vstack((Astart[0,:],Astart[2:,:],Astart[1,:]))
		A = np.hstack((A,Astart))
	
	# checks
	assert all(np.sum(Astart,axis=1) == 1)  # starter should have no overlaps
	assert all(np.sum(Astart,axis=0) == 2)  # and two pools per subject
	assert np.unique(A,axis=1).shape == A.shape  # all columns should be unique
	
	return A

def ordered_pairs(m,n):
	assert m % 2 == 0
	Abase = genpairs(m//2)
	return Abase[:,[i % Abase.shape[1] for i in range(n)]]

def grid_balanced(nlist):
	m = sum(nlist)
	n = np.product(nlist)
	q = len(nlist)
	A = np.zeros((m,n))
	i = 0
	for j1 in range(q):
		B = np.ones((1,1))
		for j2 in reversed(range(q)):
			B = np.kron(B,np.identity(nlist[j2]) if j2 == j1 else np.ones((1,nlist[j2])))
		A[i:i+nlist[j1],:] = B
		i += nlist[j1]
	return A

def partition(m,n):
	assert n/m == int(n/m)
	s = int(n/m)
	return np.kron(np.identity(m), np.ones((1,int(n/m))))

# based on R code written by Edgar
# todo - implement `read=False` case
def steiner_system(n):
	A_full = np.loadtxt("steiner_system_4_7_23.txt",skiprows=1,usecols=range(1,254))
	# select random subset
	s = np.random.choice(253,n,replace=False)
	return A_full[:,s]

def random_binary_balanced(m,n,q):
	A = np.zeros((m,n))
	for i in range(n):
		idx = np.random.choice(m,q,replace=False)
		A[idx,i] = 1
	return A

def decode(y,A,e):
	x = A.sum(0) - A.T.dot(y)
	return (x < e)

def sample_pooled_viral_load(A_i,viral_load):
	sampled_load = poisson.rvs(viral_load/A_i.sum())
	total_sampled_load = A_i.dot(sampled_load)
	return total_sampled_load

def get_n_random_samples(A,ViralLoad,LOD,fp_rate):
	idx_n = np.random.choice(ViralLoad.shape[0],A.shape[1],replace=False)
	viral_load = ViralLoad[idx_n]
	x_true = (viral_load > 0).astype(np.float)
	y = np.zeros(A.shape[0])
	for i in range(len(y)):
		total_sampled_load = sample_pooled_viral_load(A[i],viral_load)
		if (total_sampled_load > LOD) or (np.random.random() < fp_rate): # call the pool positive with probability fp_rate, regardless of viral load
			y[i] = 1
	return x_true,y,idx_n

def run_test(A,ViralLoad,InfectionTime,OnsetTime,LOD,e,fp_rate,itr=2500):
	true_pos = 0
	total_pos = 0
	T = []
	Vl = []
	Vi = []
	Vo = []
	R = []
	for __ in range(itr):
		x,y,idx_n = get_n_random_samples(A,ViralLoad,LOD,fp_rate)
		x_hat = decode(y,A,e)
		t = x_hat.sum() # all putative positives
		x_hat = x_hat*(ViralLoad[idx_n] > LOD) # samples that are putative positive and above the LOD
		true_pos += x.dot(x_hat)
		total_pos += x.sum()
		T.append(t)
		if x.sum() > 0:
			ki = np.where(x)[0]
			R.append(x[ki]*x_hat[ki])
			Vl.append(ViralLoad[idx_n][ki])
			Vi.append(InfectionTime[idx_n][ki])
			Vo.append(OnsetTime[idx_n][ki])
	T = np.array(T)
	if len(R):
		R = np.hstack(R)
		Vl = np.hstack(Vl)
		Vi = np.hstack(Vi)
		Vo = np.hstack(Vo)
	efficiency_avg = A.shape[1]/(A.shape[0] + np.average(T))
	if total_pos > 0:
		recall = true_pos/total_pos
	else:
		recall = 1
	return efficiency_avg, recall, R, Vl, Vi, Vo

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--infection-time', help='Path to time of infection for each individual')
	parser.add_argument('--onset-time', help='Path to time of symptom onset for each individual')
	parser.add_argument('--savepath', help='Path to save summary figure')
	parser.add_argument('--n-individuals', help='Number of individuals to test (n)',type=int)
	parser.add_argument('--m-pools', help='Number of pools (m)',type=int)
	parser.add_argument('--q-split', help='Number of pools per sample (q)',type=int)
	parser.add_argument('--start-time', help='Start of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	parser.add_argument('--fp-rate', help='Rate of false positive PCRs',type=float,default=0.01)
	parser.add_argument('--LOD', help='Limit of detection',type=float,default=100)
	parser.add_argument('--pool-compositions', help='Path to matrix of pool compositions, otherwise random (default)',default=None)
	parser.add_argument('--expected-pos-prob', help='1-expected faulty test rate (for setting error correction threshold)',type=float,default=0.95)
	parser.add_argument('--grid-balanced', help='List of grid dimensions, otherwise random (default)',default=None)
	parser.add_argument('--steiner', help='Number of individuals to test, otherwise random (default)',default=None)
	parser.add_argument('--ordered-pairs', help='Use ordered pairs, q-split must be 2, uses n-individuals and m-pools',action='store_true')
	parser.add_argument('--partition', help='Use a partition, q-split must be 1, uses n-individuals and m-pools',action='store_true')
	args,_ = parser.parse_known_args()
	f = open(args.viral_load_matrix)
	timepoints = f.readline()
	ViralLoad = [[float(l) for l in line.strip().split(',')] for line in f]
	f.close()
	ViralLoad = 10**np.array(ViralLoad) - 1
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
	E_avg = np.zeros(ViralLoad.shape[1])
	Recall_combined = np.zeros(ViralLoad.shape[1])
	Recall_n, Vl, Vi, Vo = [],[],[],[]
	# pooling composition matrix
	if args.pool_compositions is not None:
		print('loading pools from: %s' % args.pool_compositions)
		A = np.load(args.pool_compositions)
		m_pools = A.shape[0]
		n_individuals = A.shape[1]
		q_split = A.sum(axis=0)[0]
		assert np.all(A.sum(axis=0) == q_split)
	elif args.grid_balanced is not None:
		print('using balanced grid design: %s' % args.grid_balanced)
		nlist = [int(n) for n in args.grid_balanced.strip().split(',')]
		A = grid_balanced(nlist)
		m_pools = sum(nlist)
		n_individuals = np.product(nlist)
		q_split = len(nlist)
	elif args.steiner is not None:
		print('using steiner design: %s' % args.steiner)
		n_individuals = int(args.steiner)
		A = steiner_system(n_individuals)
		m_pools = A.shape[0]
		q_split = A.sum(axis=0)[0]
		assert np.all(A.sum(axis=0) == q_split)
	elif args.ordered_pairs:
		print('using ordered pairs design')
		m_pools = args.m_pools
		n_individuals = args.n_individuals
		q_split = args.q_split
		assert q_split == 2
		A = ordered_pairs(m_pools,n_individuals)
	elif args.partition:
		print('using partition design')
		m_pools = args.m_pools
		n_individuals = args.n_individuals
		q_split = args.q_split
		assert q_split == 1
		A = partition(m_pools,n_individuals)
	else:
		m_pools = args.m_pools
		n_individuals = args.n_individuals
		q_split = args.q_split
		A = random_binary_balanced(m_pools,n_individuals,q_split)
	# number of faulty tests to tolerate
	e = np.round(2*(1-args.expected_pos_prob)*q_split) + 1e-7
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	print('error_tol\t%d' % e)
	for t in range(args.start_time,args.end_time+1):
		ea,rec,r,vl,vi,vo = run_test(A,ViralLoad[:,t],t-InfectionTime,t-OnsetTime,args.LOD,e,args.fp_rate)
		Recall_n.append(r)
		Vl.append(vl)
		Vi.append(vi)
		Vo.append(vo)
		E_avg[t] = ea
		Recall_combined[t] = rec
		print('Time: %d; Case frequency: %.5f; Efficiency: %.2f; Sensitivity: %.4f; Fraction >LOD: %.4f' % (t,np.average(ViralLoad[:,t] > 0),ea,rec, np.average(ViralLoad[ViralLoad[:,t] > 0,t] > args.LOD)))
	np.save('%s/Eff_avg.n-%d_m-%d_q-%d.npy' % (args.savepath,n_individuals,m_pools,q_split),E_avg)
	np.save('%s/Recall_combined.n-%d_m-%d_q-%d.npy' % (args.savepath,n_individuals,m_pools,q_split),Recall_combined)
	np.save('%s/Recall_n.n-%d_m-%d_q-%d.npy' % (args.savepath,n_individuals,m_pools,q_split),Recall_n)
	np.save('%s/Vl.n-%d_m-%d_q-%d.npy' % (args.savepath,n_individuals,m_pools,q_split),Vl)
	np.save('%s/Vi.n-%d_m-%d_q-%d.npy' % (args.savepath,n_individuals,m_pools,q_split),Vi)
	np.save('%s/Vo.n-%d_m-%d_q-%d.npy' % (args.savepath,n_individuals,m_pools,q_split),Vo)

