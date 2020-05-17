import numpy as np
import argparse

# PRECOMPUTED:
# generate NP+ as a function of time relative to infection
# f = open('../../../Simulated_populations/seir_viral_loads.identification.csv')
# timepoints = f.readline()
# ViralLoad = [[float(l) for l in line.strip().split(',')] for line in f]
# f.close()
# ViralLoad = np.array(ViralLoad)
# f = open('../../../Simulated_populations/infection_times.csv')
# header = f.readline()
# InfectionTime = [float(line.strip()) for line in f]
# f.close()
# InfectionTime = np.array(InfectionTime)
# T = []
# V = []
# for i in range(len(InfectionTime)):
# 	if InfectionTime[i] > -1:
# 		idx = np.where(ViralLoad[i] > 0)[0]
# 		idx = np.arange(idx[0],ViralLoad.shape[1])
# 		T.append(idx - InfectionTime[i])
# 		V.append(ViralLoad[i,idx])
# 
# T = np.hstack(T)
# T = np.round(T,0)
# V = np.hstack(V)
# x = np.zeros(int(T.max()))
# for i in range(1,len(x)):
# 	idx = np.where(T == i)[0]
# 	x[i] = np.average(V[idx] > 2)
# 
# np.save('../../examples/study_power.NP_pos_probability.npy',x)

def exposure_probability(times,incidence,pos_prob):
	P = []
	for i in range(1,len(times)):
		inc = incidence[times[i-1]:times[i]]
		pos = pos_prob[:len(inc)][::-1]
		P.append(inc.dot(pos))
	return np.array(P)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--incidence-high', help='High incidence csv, e.g. ../../examples/study_power.incidence.high.csv')
	parser.add_argument('--incidence-low', help='Low incidence csv, e.g. ../../examples/study_power.incidence.low.csv')
	parser.add_argument('--sampling-times', help='Comma-separated list of sampling times (days), always beginning with 0, eg 0,30,60')
	parser.add_argument('--outfile', help='Output filename, e.g. ../../examples/study_power.susceptibility_table.csv')
	parser.add_argument('--NP-pos-prob', help='Path to NP+ probability relative to time of infection', default='../../examples/study_power.NP_pos_probability.npy')
	args,_ = parser.parse_known_args()
	f = open(args.incidence_high)
	header = f.readline()
	incidence_high = np.array([float(line.strip()) for line in f])
	f.close()
	f = open(args.incidence_low)
	header = f.readline()
	incidence_low = np.array([float(line.strip()) for line in f])
	f.close()
	pos_prob = np.load(args.NP_pos_prob)
	times = [int(t) for t in args.times.split(',')]
	exposure_high = exposure_probability(times,incidence_high,pos_prob)
	exposure_low = exposure_probability(times,incidence_low,pos_prob)
	IgG_neg_NP_neg_prob_high = np.product(1-exposure_high)
	IgG_neg_NP_neg_prob_low = np.product(1-exposure_low)
	S = np.linspace(0,1,100)
	IgG_pos_NP_neg_prob_high = []
	IgG_pos_NP_neg_prob_low = []
	for s in S:
		p1 = np.product(1-exposure_high*s)
		p2 = np.product(1-exposure_low*s)
		IgG_pos_NP_neg_prob_high.append(p1)
		IgG_pos_NP_neg_prob_low.append(p2)
	f = open(args.outfile,'w')
	_ = f.write(','.join(['Relative susceptibility','IgG_neg NP_pos high', 'IgG_pos NP_pos high', 'IgG_neg NP_pos low', 'IgG_pos NP_pos low']) + '\n')
	for i in range(len(S)):
		_=f.write('%f,%f,%f,%f,%f\n' % (S[i],1-IgG_neg_NP_neg_prob_high, 1-IgG_pos_NP_neg_prob_high[i], 1-IgG_neg_NP_neg_prob_low, 1-IgG_pos_NP_neg_prob_low[i]))
	f.close()

