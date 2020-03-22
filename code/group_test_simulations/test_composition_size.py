import numpy as np
import spams
import glob,os
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

def generate_random_samples(A,k,p=1):
	idx = np.random.randint(A.shape[1],size=(k,))
	x_true = np.zeros(A.shape[1])
	x_true[idx] = 1
	y = np.zeros(A.shape[0])
	for i in range(len(y)):
		z = np.random.random(A.shape[1])
		x_rand = x_true*(z < p)
		if A[i].dot(x_rand) > 0:
			y[i] = 1
	return x_true,y

def individual_pr(x_true,x_hat):
	if x_hat.sum() > 0:
		precision = x_true.dot(x_hat)/x_hat.sum()
	else:
		precision = 0
	recall = x_true.dot(x_hat)/x_true.sum()
	return precision, recall

def run_test(A,hit_rate,p,e,itr=2500):
	true_pos = 0
	total_pos = 0
	T = []
	K = []
	for __ in range(itr):
		k = (np.random.random(A.shape[1]) < hit_rate).sum()
		x,y = generate_random_samples(A,k,p)
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
	T = np.array(T)
	efficiency_avg = A.shape[1]/(A.shape[0] + np.average(T))
	efficiency_std = np.std(A.shape[1]/(A.shape[0] + T))
	return efficiency_avg, efficiency_std, true_pos/total_pos

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--center-frequency', help='Case frequency used to center the number of tests',type=float)
	parser.add_argument('--pos-prob', help='Probability of a positive result from a diluted case',type=float)
	parser.add_argument('--savepath', help='Path to save summary figure')
	parser.add_argument('--group-size-min', help='Minimum individuals for group testing range',type=int,default=50)
	parser.add_argument('--group-size-max', help='Maximum individuals for group testing range',type=int,default=400)
	parser.add_argument('--group-size-stride', help='Step size of group testing range',type=int,default=50)
	parser.add_argument('--max-dilution', help='Maximum number of samples per composition',type=int,default=100)
	parser.add_argument('--max-sample-split', help='Maximum number of compositions a single sample can be split into',type=int,default=8)
	args,_ = parser.parse_known_args()
	hr = np.logspace(np.log10(args.center_frequency/8),np.log10(args.center_frequency*8),15)
	N = np.arange(args.group_size_min,args.group_size_max,args.group_size_stride)
	E_avg = np.zeros((len(hr),len(N)))
	E_std = np.zeros((len(hr),len(N)))
	R = np.zeros((len(hr),len(N)))
	M = np.zeros(len(N))
	for j,n in enumerate(N):
		m0 = int(np.log2(n)*max((n*args.center_frequency),1)/args.pos_prob**3)
		m_range = np.unique(np.linspace(max(1,m0/4),m0*1.5,10,dtype=np.int))
		e_avg = np.zeros((len(hr),len(m_range)))
		e_std = np.zeros((len(hr),len(m_range)))
		r = np.zeros((len(hr),len(m_range)))
		for i,hit_rate_actual in enumerate(hr):
			for l,m in enumerate(m_range):
				q = int(min(args.max_dilution/n*m,args.max_sample_split,m/2)) # could also consider np.log(n) upper bound here
				A = random_binary_balanced(m,n,q)
				e = max(1.5*(1-args.pos_prob)*np.average(A)*A.shape[0], 1.1)
				result = run_test(A,hit_rate_actual,args.pos_prob,e)
				e_avg[i,l] = result[0]
				e_std[i,l] = result[1]
				r[i,l] = result[2]
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
	plt.savefig('%s/efficiency.center-%.4f_pos-prob-%.2f.png' % (args.savepath,args.center_frequency,args.pos_prob))
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
	plt.savefig('%s/recall.center-%.4f_pos-prob-%.2f.png' % (args.savepath,args.center_frequency,args.pos_prob))
	plt.close()
	np.save('%s/Eff_avg.center-%.4f_pos-prob-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob),E_avg)
	np.save('%s/Eff_std.center-%.4f_pos-prob-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob),E_std)
	np.save('%s/Recall.center-%.4f_pos-prob-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob),R)
	np.save('%s/M.center-%.4f_pos-prob-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob),M)
	np.save('%s/N.center-%.4f_pos-prob-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob),N)
	np.save('%s/X.center-%.4f_pos-prob-%.2f.npy' % (args.savepath,args.center_frequency,args.pos_prob),hr)



