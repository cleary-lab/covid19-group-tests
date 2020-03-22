import numpy as np
import spams
import glob,os
from scipy.stats import truncnorm
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns
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

def generate_random_samples(A,k,p=1,s=0):
	if s > 0:
		# 1-faulty test probabilities per individual drawn from a truncated normal distribution
		p1 = truncnorm.rvs((0-p)/s,(1-p)/s,loc=p,scale=s,size=A.shape[1])
	else:
		# constant faulty test probability
		p1 = p
	# randomly assign k positive samples to n individuals
	idx = np.random.randint(A.shape[1],size=(k,))
	x_true = np.zeros(A.shape[1])
	x_true[idx] = 1
	y = np.zeros(A.shape[0])
	for i in range(len(y)):
		z = np.random.random(A.shape[1])
		# randomly generate faults for each individual in this pool
		x_rand = x_true*(z < p1)
		if A[i].dot(x_rand) > 0:
			# pool returned positive result
			y[i] = 1
	return x_true,y,p1

def run_test(A,hit_rate,p,e,s=0,itr=2500):
	true_pos = 0
	total_pos = 0
	T = []
	K = []
	P = []
	R = []
	for __ in range(itr):
		# random number of positive samples
		k = (np.random.random(A.shape[1]) < hit_rate).sum()
		x,y,p1 = generate_random_samples(A,k,p=p,s=s)
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
		K.append(k)
		if (k > 0) and (s > 0):
			ki = np.where(x)[0]
			R.append(x[ki]*x_hat3[ki])
			P.append(p1[ki])
	T = np.array(T)
	if s > 0:
		R = np.hstack(R)
		P = np.hstack(P)
	efficiency_avg = A.shape[1]/(A.shape[0] + np.average(T))
	efficiency_std = np.std(A.shape[1]/(A.shape[0] + T))
	return efficiency_avg, efficiency_std, true_pos/total_pos, R, P

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--center-frequency', help='Case frequency used to center the number of tests',type=float)
	parser.add_argument('--pos-prob', help='Probability mean (intruncated normal) of a positive result from a diluted case',type=float)
	parser.add_argument('--pos-prob-std', help='Stdev (in truncated normal) of positive probability',type=float,default=0.2)
	parser.add_argument('--savepath', help='Path to save summary figure')
	parser.add_argument('--group-size-min', help='Minimum individuals (n) for group testing range',type=int,default=50)
	parser.add_argument('--group-size-max', help='Maximum individuals (n) for group testing range',type=int,default=400)
	parser.add_argument('--group-size-stride', help='Step size of group testing range',type=int,default=50)
	parser.add_argument('--max-dilution', help='Maximum number of samples per composition',type=int,default=100)
	parser.add_argument('--max-sample-split', help='Maximum number of compositions a single sample can be split into',type=int,default=8)
	args,_ = parser.parse_known_args()
	hr = np.logspace(np.log10(args.center_frequency/8),np.log10(args.center_frequency*8),15)
	# test a range for values of n (total number of individuals to be split into pools)
	N = np.arange(args.group_size_min,args.group_size_max,args.group_size_stride)
	E_avg = np.zeros((len(hr),len(N)))
	E_std = np.zeros((len(hr),len(N)))
	R = np.zeros((len(hr),len(N)))
	M = np.zeros(len(N))
	p1 = truncnorm.rvs((0-args.pos_prob)/args.pos_prob_std,(1-args.pos_prob)/args.pos_prob_std,loc=args.pos_prob,scale=args.pos_prob_std,size=1000)
	for j,n in enumerate(N):
		m0 = np.round(np.average(np.log2(n)*max((n*args.center_frequency),1)/p1**3)).astype(np.int)
		# test a range for m (number of pools)
		m_range = np.unique(np.linspace(max(1,m0/3),m0*1.5,10,dtype=np.int))
		e_avg = np.zeros((len(hr),len(m_range)))
		e_std = np.zeros((len(hr),len(m_range)))
		r = np.zeros((len(hr),len(m_range)))
		for i,hit_rate_actual in enumerate(hr):
			for l,m in enumerate(m_range):
				if m < n/4:
					# number of pools into which each sample is split
					q = int(min(args.max_dilution/n*m,args.max_sample_split,m/2,np.round(np.log(n))))
					# pooling composition matrix
					A = random_binary_balanced(m,n,q)
					# number of faulty tests to tolerate
					e = np.round(2*(1-np.average(p1))*np.average(A)*A.shape[0]) + 1e-7
					ea,es,rec,rr,pp = run_test(A,hit_rate_actual,args.pos_prob,e,s=args.pos_prob_std)
					e_avg[i,l] = ea
					e_std[i,l] = es
					r[i,l] = rec
				else:
					e_avg[i,l] = 1
					e_std[i,l] = 0
					r[i,l] = 1
		best_idx = np.argmax(np.average(e_avg*r,axis=0))
		E_avg[:,j] = e_avg[:,best_idx]
		E_std[:,j] = e_std[:,best_idx]
		R[:,j] = r[:,best_idx]
		M[j] = m_range[best_idx]
		print(n,np.average(E_avg[:,j]),np.average(R[:,j]),M[j])
	for n in N:
		j = np.where(N == n)[0][0]
		label = 'n=%d, m=%d, avg_recall=%.3f' % (n,M[j],np.average(R[:,j]))
		_=plt.plot(np.log(hr), E_avg[:,j], label=label)
	_=plt.xticks(np.log([hr.min(),hr.min()*5,hr.max()/5,hr.max()]),[np.around(hr.min(),4),np.around(hr.min()*5,4),np.around(hr.max()/5,4),np.around(hr.max(),4)])
	_=plt.legend(loc='best',prop={'size': 8})
	_=plt.xlabel('Frequency of cases')
	_=plt.ylabel('Efficiency (# individuals / # tests)')
	_=plt.title('Tests designed assuming average case frequency of %.3f\nProbability of positive result from a case: %.2f' % (args.center_frequency,args.pos_prob))
	plt.savefig('%s/efficiency.center-%.4f_pos-prob-%.2f_prob-std-%.2f.png' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std))
	plt.close()
	for n in N:
		j = np.where(N == n)[0][0]
		label = 'n=%d, m=%d, avg_efficiency=%.3f' % (n,M[j],np.average(E_avg[:,j]))
		_=plt.plot(np.log(hr), R[:,j], label=label)
	_=plt.xticks(np.log([hr.min(),hr.min()*5,hr.max()/5,hr.max()]),[np.around(hr.min(),4),np.around(hr.min()*5,4),np.around(hr.max()/5,4),np.around(hr.max(),4)])
	_=plt.legend(loc='best',prop={'size': 8})
	_=plt.xlabel('Frequency of cases')
	_=plt.ylabel('Recall')
	_=plt.title('Tests designed assuming average case frequency of %.3f\nProbability of positive result from a case: %.2f' % (args.center_frequency,args.pos_prob))
	plt.savefig('%s/recall.center-%.4f_pos-prob-%.2f_prob-std-%.2f.png' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std))
	plt.close()
	np.save('%s/Eff_avg.center-%.4f_pos-prob-%.2f_prob-std-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std),E_avg)
	np.save('%s/Eff_std.center-%.4f_pos-prob-%.2f_prob-std-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std),E_std)
	np.save('%s/Recall.center-%.4f_pos-prob-%.2f_prob-std-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std),R)
	np.save('%s/M.center-%.4f_pos-prob-%.2f_prob-std-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std),M)
	np.save('%s/N.center-%.4f_pos-prob-%.2f_prob-std-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std),N)
	np.save('%s/X.center-%.4f_pos-prob-%.2f_prob-std-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob,args.pos_prob_std),hr)



