import numpy as np
from scipy.special import comb
from scipy.stats import poisson
from sklearn.neighbors import KernelDensity
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')

def empirical_convolutions(ViralLoad,K,b=0.1,samples=10000,thresh=0):
	Kernels = []
	nz = np.where(ViralLoad.flatten() > thresh)[0]
	x = ViralLoad.flatten()[nz]
	x = 10**x
	for k in range(1,K+1):
		z = [np.random.choice(x,k).sum() for _ in range(samples)]
		z = np.log10(z)
		kernel = KernelDensity(bandwidth=b).fit(z[:,np.newaxis])
		Kernels.append(kernel)
	return Kernels

def total_density(m,p,conditional_k):
	x = 0
	for k in range(m+1):
		x += comb(m,k)*p**k*(1-p)**(m-k)*conditional_k[k]
	return x

def cond_prob(p,m,k,td,conditional_k):
	cp = comb(m,k)*p**k*(1-p)**(m-k)*conditional_k[k]/td
	return cp

def run_em(y,b,n,p0,Kernels,max_itr=100,min_itr=10,tol=1e-3,n_vals=50):
	# calculate all the conditional densities for the observed values (y) up front
	CK = []
	for i in range(b):
		if y[i] > 0:
			conditional_k = [0]+[np.exp(Kernels[k].score_samples([[np.log10(y[i]*n)]])[0]) for k in range(n)] # conditional density of y*m, given k
		else:
			conditional_k = [1]+[np.exp(Kernels[k].score_samples([[-1]])[0]) for k in range(n)] # conditional density of 0.1
		CK.append(conditional_k)
	for itr in range(max_itr):
		p = 0
		for i in range(b):
			conditional_k = CK[i]
			td = total_density(n,p0,conditional_k)
			for k in range(1,n+1):
				p += k*cond_prob(p0,n,k,td,conditional_k)
		p = p/(n*b)
		if (itr > min_itr) and (abs(p-p0)/p0 < tol):
			break
		p0 = p
	return p

def sample_pooled_viral_load(freq,nonzero_viral_loads,n_samples,LOD=1):
	z = np.random.random(n_samples)
	n = (z < freq).sum()
	if n > 0:
		probs = np.array([len(vl) for vl in nonzero_viral_loads]) # sample individuals with longer infection periods more often
		probs = probs/probs.sum()
		idx = np.random.choice(len(nonzero_viral_loads),n,replace=False,p=probs)
		vl = np.array([np.random.choice(nonzero_viral_loads[i],1) for i in idx])
		vl = 10**vl - 1
		sampled_load = poisson.rvs(vl/n_samples)
		if n > 1:
			total_sampled_load = sampled_load.sum()
		else:
			total_sampled_load = sampled_load
		if total_sampled_load < LOD:
			total_sampled_load = 0
	else:
		total_sampled_load = 0
	sampled_freq = n/n_samples
	return sampled_freq,total_sampled_load


f = open('viral_loads.txt')
ViralLoad = [[float(l) for l in line.strip().split()] for line in f]
f.close()
ViralLoad = np.array(ViralLoad)

N=24*384 # total individuals
b=24 # number of batches
n=384 # size of each batch

Kernels = empirical_convolutions(ViralLoad,n)

nonzero_viral_loads = [vl[vl > 0] for vl in ViralLoad if vl.max()>1]

frequency_range = 10**np.linspace(-4,-1,25)

P = []
for freq in frequency_range:
	result = [sample_pooled_viral_load(freq,nonzero_viral_loads,n) for _ in range(n)]
	Y = [r[1] for r in result]
	p1 = np.average([r[0] for r in result])
	if max(Y) > 0:
		p_hat = run_em(Y,b,n,0.005,Kernels)
	else:
		p_hat = 0
	P.append((p1,p_hat))
	print(freq,p1,p_hat)


P = np.log10(P)
P[P < np.log10(1/N)] = np.log10(1/N)
_=plt.plot(np.log10(frequency_range),np.log10(frequency_range),label='entire population; N=%d' % ViralLoad.shape[0])
_=plt.plot(np.log10(frequency_range),P[:,0],label='sampled population; N=%d,b=%d,n=%d' % (N,b,n))
_=plt.plot(np.log10(frequency_range),P[:,1],label='estimated from b=%d measurements' % b)
_=plt.legend(loc='best',prop={'size': 8})
_=plt.xlabel('log( true frequency of positive samples )',fontsize=10)
_=plt.ylabel('log( estimated frequency of positive samples )',fontsize=9)
_=plt.tight_layout()
plt.savefig('estimated_prevalence.b-%d_n-%d.png' % (b,n))
plt.close()

