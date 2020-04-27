import numpy as np
import argparse
from scipy.special import comb
from scipy import integrate
import joblib
import sys

def total_density(m,p,conditional_k):
	x = 0
	for k in range(m+1):
		x += comb(m,k)*p**k*(1-p)**(m-k)*conditional_k[k]
	return x

def cond_prob(p,m,k,td,conditional_k):
	cp = comb(m,k)*p**k*(1-p)**(m-k)*conditional_k[k]/td
	return cp

def undetected_prob(Kernels,n,fp_rate,Ct_thresh=40):
	# conditional probability of Ct value > Ct_thresh
	lower = Ct_thresh
	upper = np.inf
	ck0 = [1-fp_rate]
	for k in range(n):
		# integrate over Ct values > Ct_thresh
		f = lambda x: np.exp(Kernels[k].score_samples([[x]])[0])
		p = integrate.quad(f,lower,upper)[0]
		ck0.append(p)
	return ck0

def run_em(y,b,n,p0,Kernels,fp_rate,max_itr=100,min_itr=10,tol=1e-3,n_vals=50):
	# calculate all the conditional densities for the observed values (y) up front
	CK = []
	ck0 = undetected_prob(Kernels,n,fp_rate)
	for i in range(b):
		if y[i] > -1: # -1 indicates undetected
			conditional_k = [np.exp(Kernels[0].score_samples([[y[i] - np.log2(10)*np.log10(n)]])[0])*fp_rate]+[np.exp(Kernels[k].score_samples([[y[i] - np.log2(10)*np.log10(n)]])[0]) for k in range(n)] #conditional density of y*n, given k
		else:
			conditional_k = ck0
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

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--batch-size', help='Number of samples in each batch',type=int)
	parser.add_argument('--num-batches', help='Number of batches',type=int)
	parser.add_argument('--observations', help='Path to file with viral load observed in each batch')
	parser.add_argument('--kde-path', help='Path to KDEs fit for each value of k, including prefix (ie, seir_kde)')
	parser.add_argument('--p0', help='Initial frequency estimate',type=float,default=0.005)
	parser.add_argument('--fp-rate', help='Estimated false positive PCR rate',type=float,default=0.002)
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	f = open(args.observations)
	Y = [float(line.strip()) for line in f]
	f.close()
	if len(Y) != args.num_batches:
		print('Number of observations in %s does not equal batch size' % args.observations)
		sys.exit()
	Kernels = [joblib.load('%s.k-%d.pkl' % (args.kde_path,k)) for k in range(1,args.batch_size+1)]
	p_est = run_em(Y,args.num_batches,args.batch_size,args.p0,Kernels,args.fp_rate)
	print('Estimated prevalence: %.4f%%' % (p_est*100))
