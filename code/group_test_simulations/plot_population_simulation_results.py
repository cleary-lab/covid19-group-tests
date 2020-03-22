import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns

def bin_results(R,P,n=50):
	xs = np.argsort(P)
	bins = [P[xs[i]] for i in range(0,len(P),int(len(P)/n))]
	indices = np.digitize(P,bins)
	x = np.array([np.average(P[indices == i]) for i in set(indices)])
	y = np.array([np.average(R[indices == i]) for i in set(indices)])
	return x,y

f = open('../tmp.txt')
ViralLoad = [[float(l) for l in line.strip().split()] for line in f]
f.close()
ViralLoad = np.array(ViralLoad)
x = [np.average(v[v > 0] > 2) for v in ViralLoad.T]

recall = np.load('Recall_combined.n-1000_m-30.npy')
efficiency = np.load('Eff_avg.n-1000_m-30.npy')
recall_n = np.load('Recall_n.n-1000_m-30.npy',allow_pickle=True)
Vl = np.load('Vl.n-1000_m-30.npy',allow_pickle=True)
Vt = np.load('Vt.n-1000_m-30.npy',allow_pickle=True)

v1,r1 = bin_results(np.hstack(recall_n), np.hstack(Vl))
v2,r2 = bin_results(np.hstack(recall_n), np.hstack(Vt))

ax1 = plt.subplot(2,2,1)
_=plt.plot(recall,label='group test')
_=plt.plot(x,label='individual test')
_=plt.legend(loc='best',prop={'size': 8})
_=plt.xlabel('time')
_=plt.ylabel('Average recall')

ax2 = plt.subplot(2,2,3)
_=plt.plot(efficiency)
_=plt.xlabel('Time')
_=plt.ylabel('Average efficiency')

ax3 = plt.subplot(2,2,2)
_=plt.plot(np.log10(v1),r1)
_=plt.xlabel('log10( Viral load )')
_=plt.ylabel('Average recall')

ax4 = plt.subplot(2,2,4)
_=plt.scatter(v2-5,r2)
_=plt.xlabel('Time since symptom onset')
_=plt.ylabel('Average recall')

_=plt.tight_layout()
plt.savefig('../tmp.png')
plt.close()


