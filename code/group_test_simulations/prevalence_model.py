import numpy as np
from scipy.special import comb
from scipy.stats import poisson
from sklearn.neighbors import KernelDensity
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')

def empirical_distribution(ViralLoad,b=0.1):
	nz = np.where(ViralLoad.flatten() > 0)[0]
	x = ViralLoad.flatten()[nz]
	kernel = KernelDensity(bandwidth=b).fit(x[:,np.newaxis])
	return kernel

def convolve_empirical_k(density,k):
	z = np.copy(density)
	for i in range(1,k):
		z = np.convolve(density,z)
	return z

def total_density(m,p,convolve_k):
	x = 0
	for k in range(m+1):
		x += comb(m,k)*p**k*(1-p)**(m-k)*convolve_k[k]
	return x

def cond_prob(p,m,k,td,convolve_k):
	cp = comb(m,k)*p**k*(1-p)**(m-k)*convolve_k[k]/td
	return cp

def run_em(y,n,m,p0,kernel,max_itr=100,min_itr=10,tol=1e-3):
	CK = []
	for i in range(n):
		if y[i] > 0:
			density_vals = np.log10(np.linspace(1,y[i]*m,100))[:,np.newaxis]
			density = np.exp(kernel.score_samples(density_vals))
			convolve_k = [0]+[convolve_empirical_k(density,k)[99] for k in range(1,m+1)] # conditional density of y*m, given k
		else:
			convolve_k = [1]+[np.exp(kernel.score_samples([[-1]])[0])**k for k in range(1,m+1)]
		CK.append(convolve_k)
	for itr in range(max_itr):
		p = 0
		for i in range(n):
			convolve_k = CK[i]
			td = total_density(m,p0,convolve_k)
			for k in range(1,m+1):
				p += k*cond_prob(p0,m,k,td,convolve_k)
		p = p/(n*m)
		if (itr > min_itr) and (abs(p-p0)/p0 < tol):
			break
		p0 = p
	return p

def sample_pooled_viral_load(freq,nonzero_viral_loads,n_samples,LOD=1):
	z = np.random.random(n_samples)
	n = (z < freq).sum()
	if n > 0:
		idx = np.random.choice(len(nonzero_viral_loads),n,replace=False)
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
n=24 # number of batches
m=int(N/n) # size of each batch

kernel = empirical_distribution(ViralLoad,b=0.1)
nonzero_viral_loads = [vl[vl > 0] for vl in ViralLoad if vl.max()>1]

frequency_range = 10**np.linspace(-4,-1,25)

P = []
for freq in frequency_range:
	result = [sample_pooled_viral_load(freq,nonzero_viral_loads,m) for _ in range(n)]
	Y = [r[1] for r in result]
	p1 = np.average([r[0] for r in result])
	if max(Y) > 0:
		p_hat = run_em(Y,n,m,0.005,kernel)
	else:
		p_hat = 0
	P.append((p1,p_hat))
	print(freq,p1,p_hat)


P = np.log10(P)
P[P < np.log10(1/(n*m))] = np.log10(1/(n*m))
_=plt.plot(np.log10(frequency_range),np.log10(frequency_range),label='entire population; N=%d' % ViralLoad.shape[0])
_=plt.plot(np.log10(frequency_range),P[:,0],label='sampled population; N=%d,n=%d,m=%d' % (N,n,m))
_=plt.plot(np.log10(frequency_range),P[:,1],label='estimated from n=%d measurements' % n)
_=plt.legend(loc='best',prop={'size': 8})
_=plt.xlabel('log( true frequency of positive samples )',fontsize=10)
_=plt.ylabel('log( estimated frequency of positive samples )',fontsize=9)
_=plt.tight_layout()
plt.savefig('estimated_prevalence.png')
plt.close()

