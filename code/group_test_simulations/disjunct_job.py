import numpy as np
import sys

def random_binary_balanced(m,n,q):
	A = np.zeros((m,n))
	for i in range(n):
		idx = np.random.choice(m,q,replace=False)
		A[idx,i] = 1
	return A

def average_disjunctness(A,d,t=20000):
	z = 0
	for _ in range(t):
		s = np.random.choice(A.shape[1],d,replace=False)
		u = np.setdiff1d(np.arange(A.shape[1]),s)
		x = A[:,u].T.dot(A[:,s].max(1)) - A[:,u].sum(0)
		if x.max() < 0:
			z += 1
	return z/t


m,n,q,d,t,job_id = sys.argv[1:]
m = int(m)
n = int(n)
q = int(q)
d = int(d)
t = int(t)

best = 0
for _ in range(t):
	a = random_binary_balanced(m,n,q)
	disj = average_disjunctness(a,d)
	if disj > best:
		best = disj
		A = a

np.save('/broad/clearylab/GroupTesting/search_m-%d_n-%d_q-%d_d-%d/%s.d.npy' % (m,n,q,d,job_id),best)
np.save('/broad/clearylab/GroupTesting/search_m-%d_n-%d_q-%d_d-%d/%s.A.npy' % (m,n,q,d,job_id),A)