import numpy as np
import sys
import spams
import glob,os


def random_binary_balanced(m,n,q):
	A = np.zeros((m,n))
	for i in range(n):
		idx = np.random.choice(m,q,replace=False)
		A[idx,i] = 1
	return A

def decode(y,A,e):
	x = A.sum(0) - A.T.dot(y)
	return (x < e)

THREADS=1
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



m,n,hit_rate_actual,pos_prob,t,job_id = sys.argv[1:]
m = int(m)
n = int(n)
hit_rate_actual = float(hit_rate_actual)
pos_prob = float(pos_prob)
t = int(t)
max_dilution = 100
max_sample_split = 8

best_x = 0
for _ in range(t):
	q = int(min(max_dilution/n*m,max_sample_split,m/2))
	A = random_binary_balanced(m,n,q)
	e = max(1.5*(1-pos_prob)*np.average(A)*A.shape[0], 1.1)
	result = run_test(A,hit_rate_actual,pos_prob,e)
	x = result[0]*result[2]
	if x > best_x:
		best_x = x
		best_r = (result[0],result[2])
		best_A = A


np.save('/broad/clearylab/GroupTesting/search_m-%d_n-%d_hr-%.4f_p-%.2f/%s.r.npy' % (m,n,hit_rate_actual,pos_prob,job_id),best_r)
np.save('/broad/clearylab/GroupTesting/search_m-%d_n-%d_hr-%.4f_p-%.2f/%s.A.npy' % (m,n,hit_rate_actual,pos_prob,job_id),A)