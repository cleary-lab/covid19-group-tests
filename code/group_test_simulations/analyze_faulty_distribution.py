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
		p1 = truncnorm.rvs((0-p)/s,(1-p)/s,loc=p,scale=s,size=A.shape[1])
	else:
		p1 = p
	idx = np.random.randint(A.shape[1],size=(k,))
	x_true = np.zeros(A.shape[1])
	x_true[idx] = 1
	y = np.zeros(A.shape[0])
	for i in range(len(y)):
		z = np.random.random(A.shape[1])
		x_rand = x_true*(z < p1)
		if A[i].dot(x_rand) > 0:
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

def bin_results(R,P,n=50):
	xs = np.argsort(P)
	bins = [P[xs[i]] for i in range(0,len(P),int(len(P)/n))] + [1]
	indices = np.digitize(P,bins)
	x = np.array([np.average(P[indices == i]) for i in set(indices)])
	y = np.array([np.average(R[indices == i]) for i in set(indices)])
	return x,y


m,n,q,p,s= 35,200,5,0.9,0.2
p1 = np.average(truncnorm.rvs((0-p)/s,(1-p)/s,loc=p,scale=s,size=1000))
A = random_binary_balanced(m,n,q)
e = np.round(2.55*(1-p1)*np.average(A)*A.shape[0]) + 1e-7
ea,es,rec,R,P = run_test(A,0.0031,p,e,s=s,itr=15000)
x,y = bin_results(R,P)
print(A.shape,np.average(A.sum(1)),q,e)
print(ea,rec)

p1 = truncnorm.rvs((0-p)/s,(1-p)/s,loc=p,scale=s,size=1000)
n=750
m = np.round(np.average(np.log2(n)*max((n*0.005),1)/p1**3)).astype(np.int)
for m in range(30,121,10):
	q = int(min(100/n*m,8,m/2,np.round(np.log(n))))
	A = random_binary_balanced(m,n,q)
	e = np.round(2*(1-np.average(p1))*np.average(A)*A.shape[0]) + 1e-7
	ea,es,rec,R,P = run_test(A,0.0006,p,e,s=s,itr=15000)
	x,y = bin_results(R,P)
	print(A.shape,np.average(A.sum(1)),q,e)
	print(ea,rec)
	print('')



for d in np.linspace(1.5,10,15):
	z = np.exp(-d**2*(1-p1)*np.average(A)*A.shape[0]/(d+2))
	print(d,z)

ax1=plt.subplot(2,1,1)
_=plt.hist(P,bins=25,density=True, histtype='stepfilled', alpha=0.2)
_=plt.ylabel('density')
ax2=plt.subplot(2,1,2)
_=plt.plot(x,y)
_=plt.ylabel('Average Recall')
_=plt.xlabel('1 - faulty test probability')
_=ax2.set_xlim(ax1.get_xlim())
_=plt.savefig('tmp.png')
plt.close()

