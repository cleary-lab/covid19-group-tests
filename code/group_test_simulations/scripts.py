import numpy as np
import spams
import glob,os

def random_binary(m,n,q):
	A = np.random.random(m*n)
	return (A < q).reshape((m,n)).astype(np.float)

def random_binary_balanced(m,n,q):
	A = np.zeros((m,n))
	for i in range(n):
		idx = np.random.choice(m,q,replace=False)
		A[idx,i] = 1
	return A

def disjunctness(A):
	w = A.sum(0).min()
	z = (A.T.dot(A) - np.eye(A.shape[1])).max()
	return (w-1)/z

def average_disjunctness(A,d,t=2000):
	z = 0
	for _ in range(t):
		s = np.random.choice(A.shape[1],d,replace=False)
		u = np.setdiff1d(np.arange(A.shape[1]),s)
		x = A[:,u].T.dot(A[:,s].max(1)) - A[:,u].sum(0)
		if x.max() < 0:
			z += 1
	return z/t

def decode(y,A,e):
	# a naive decoder
	# based on the set difference: {tests containing sample i} - {tests with positive hits}
	# if this difference is less than e, return i as a possible hit
	x = A.sum(0) - A.T.dot(y)
	return (x < e)

THREADS=40
def sparse_decode(Y,D,lda):
	Ynorm = np.linalg.norm(Y)**2/Y.shape[1]
	W = spams.lasso(np.asfortranarray(Y),np.asfortranarray(D),lambda1=lda*Ynorm,mode=1,numThreads=THREADS,pos=True)
	W = np.asarray(W.todense())
	return W

def pcr_pos_prob():
	# x: log( copy number )
	# b: (cycle number)
	# c: log( copy number ) for 50% probability of detection
	f = 1/(1+np.exp(-b*(x-c)))

def generate_random_samples(A,k,p=1,s=0):
	if s > 0:
		p = np.random.randn(A.shape[1])*s + p
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

def make_our_rs_mat(N=500,m=5,q=7):
	deg = int(log( N, q ))  # make sure we have enough columns to cover N
	# this is the ordering that will order polynomials like "0, 1, 2, .. ,15, x, x+1, x+2, .. ,x+15, 2x, 2x + 1, ...."
	C1,image,polys = create_rs_code( q, deg, concatenated=True, orderoption="modoutconstant", orderfieldelts="lex")
	return C1[ :m*q, :N ]

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

def do_q(q,dil):
	m0 = q*2
	for i in range(1,20):
		print(m0*i,dil*m0*i/q)


FP = glob.glob(os.path.join('../search_m-18_n-450_q-4_d-1/*.d.npy'))
D = [np.load(fp) for fp in FP]
#best_idx = np.argmax(np.product(D,axis=1))
best_idx = np.argmax(D)
D[best_idx]
A = np.load(FP[best_idx].replace('.d.','.A.'))
average_disjunctness(A,1)
average_disjunctness(A,2)
average_disjunctness(A,3)

hit_rate = 0.0002

n=500
p = 0.90
m = int(np.log2(n)*max((n*hit_rate),1)/p**3)
q = int(min(100/n*m,8,m/2))
A = random_binary_balanced(m,n,q)
e = max(1.5*(1-p)*np.average(A)*A.shape[0], 1.1)
print(A.shape,np.average(A.sum(1)),q,e)

hit_rate = 0.0002

P1 = []
R1 = []
P2 = []
R2 = []
P3 = []
R3 = []
T = []
K = []
true_pos = 0
total_pos = 0
for __ in range(2500):
	k = (np.random.random(A.shape[1]) < hit_rate).sum()
	x,y = generate_random_samples(A,k,p)
	x_hat1 = decode(y,A,e)
	x_hat2 = sparse_decode(y[:,np.newaxis],A,0.01)
	x_hat2 = (x_hat2[:,0] > 0.01)
	x_hat3 = ((x_hat1 + x_hat2) > 0)
	true_pos += x.dot(x_hat3)
	total_pos += x.sum()
	if k > 0:
		p1,r1 = individual_pr(x,x_hat1)
		p2,r2 = individual_pr(x,x_hat2)
		p3,r3 = individual_pr(x,x_hat3)
		P1.append(p1)
		R1.append(r1)
		P2.append(p2)
		R2.append(r2)
		P3.append(p3)
		R3.append(r3)
	T.append(x_hat3.sum())
	K.append(k)

np.average(P1),np.average(R1)
np.average(P2),np.average(R2)
np.average(P3),np.average(R3)
np.average(K)
true_pos/total_pos
A.shape[1]/(A.shape[0] + np.average(T))


.77, 5.73	0.77,5.96
.97, 2.98	0.97,2.92


python test_composition_size.py --center-frequency 0.03 --pos-prob 0.9 --group-size-min 25 --group-size-max 155 --group-size-stride 25 --savepath ../Results
python test_composition_size.py --center-frequency 0.03 --pos-prob 0.8 --group-size-min 25 --group-size-max 155 --group-size-stride 25 --savepath ../Results
python test_composition_size.py --center-frequency 0.03 --pos-prob 0.7 --group-size-min 25 --group-size-max 155 --group-size-stride 25 --savepath ../Results

python test_composition_size.py --center-frequency 0.01 --pos-prob 0.9 --group-size-min 30 --group-size-max 245 --group-size-stride 30 --savepath ../Results
python test_composition_size.py --center-frequency 0.01 --pos-prob 0.8 --group-size-min 30 --group-size-max 245 --group-size-stride 30 --savepath ../Results
python test_composition_size.py --center-frequency 0.01 --pos-prob 0.7 --group-size-min 30 --group-size-max 245 --group-size-stride 30 --savepath ../Results

python test_composition_size.py --center-frequency 0.005 --pos-prob 0.9 --group-size-min 50 --group-size-max 355 --group-size-stride 50 --savepath ../Results
python test_composition_size.py --center-frequency 0.005 --pos-prob 0.8 --group-size-min 50 --group-size-max 355 --group-size-stride 50 --savepath ../Results
python test_composition_size.py --center-frequency 0.005 --pos-prob 0.7 --group-size-min 50 --group-size-max 355 --group-size-stride 50 --savepath ../Results

python test_composition_size.py --center-frequency 0.001 --pos-prob 0.9 --group-size-min 100 --group-size-max 455 --group-size-stride 50 --savepath ../Results
python test_composition_size.py --center-frequency 0.001 --pos-prob 0.8 --group-size-min 100 --group-size-max 455 --group-size-stride 50 --savepath ../Results
python test_composition_size.py --center-frequency 0.001 --pos-prob 0.7 --group-size-min 100 --group-size-max 455 --group-size-stride 50 --savepath ../Results


python test_composition_size.py --center-frequency 0.03 --pos-prob 0.95 --group-size-min 25 --group-size-max 155 --group-size-stride 25 --savepath ../Results
python test_composition_size.py --center-frequency 0.01 --pos-prob 0.95 --group-size-min 30 --group-size-max 245 --group-size-stride 30 --savepath ../Results
python test_composition_size.py --center-frequency 0.005 --pos-prob 0.95 --group-size-min 50 --group-size-max 355 --group-size-stride 50 --savepath ../Results
python test_composition_size.py --center-frequency 0.001 --pos-prob 0.95 --group-size-min 100 --group-size-max 455 --group-size-stride 50 --savepath ../Results



python test_composition_size.py --center-frequency 0.005 --pos-prob 0.9 --pos-prob-std 0.2 --group-size-min 300 --group-size-max 901 --group-size-stride 100 --savepath ../Results_variableP
python test_composition_size.py --center-frequency 0.005 --pos-prob 0.9 --pos-prob-std 0.2 --group-size-min 50 --group-size-max 251 --group-size-stride 50 --savepath ../Results_variableP2




