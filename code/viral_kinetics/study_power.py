import numpy as np
import argparse

# generate 2 incidence curves using simulate_data.R
# curve 1
seir_pars <- c("R0"=1.5,"gamma"=1/7,"sigma"=1/6.4,"I0"= population_n*0.005, "recovered0"=population_n*.1)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
pdf('../../examples/study_power.incidence.high.pdf')
epidemic_process$plot
dev.off()
# 	incidence = viral_loads <- as_tibble(epidemic_process$incidence)
# 	write_csv(incidence,"../../examples/study_power.incidence.high.csv")
# 
# curve 2
seir_pars <- c("R0"=1.1,"gamma"=1/7,"sigma"=1/6.4,"I0"= population_n*0.005, "recovered0"=population_n*.2)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
pdf('../../examples/study_power.incidence.low.pdf')
epidemic_process$plot
dev.off()
# 	incidence = viral_loads <- as_tibble(epidemic_process$incidence)
# 	write_csv(incidence,"../../examples/study_power.incidence.low.csv")

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


f = open('../../examples/study_power.incidence.high.csv')
header = f.readline()
incidence_high = np.array([float(line.strip()) for line in f])
f.close()
f = open('../../examples/study_power.incidence.low.csv')
header = f.readline()
incidence_low = np.array([float(line.strip()) for line in f])
f.close()
pos_prob = np.load('../../examples/study_power.NP_pos_probability.npy')

times = [0,30,60,90,120,150,180]
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

f = open('../../examples/study_power.susceptibility_table.csv','w')
_ = f.write(','.join(['Relative susceptibility','IgG_neg NP_pos high', 'IgG_pos NP_pos high', 'IgG_neg NP_pos low', 'IgG_pos NP_pos low']) + '\n')
for i in range(len(S)):
	_=f.write('%f,%f,%f,%f,%f\n' % (S[i],1-IgG_neg_NP_neg_prob_high, 1-IgG_pos_NP_neg_prob_high[i], 1-IgG_neg_NP_neg_prob_low, 1-IgG_pos_NP_neg_prob_low[i]))

f.close()


# back in R
pow_func = function(pars){
p1 = power.fisher.test(pars$IgG_pos.NP_pos.high, pars$IgG_neg.NP_pos.high, n1, n2, alpha=0.05, nsim=100, alternative="two.sided")
p2 = power.fisher.test(pars$IgG_pos.NP_pos.low, pars$IgG_neg.NP_pos.low, n1, n2, alpha=0.05, nsim=100, alternative="two.sided")
return(list(p1,p2))
}

S = read.table('../../examples/study_power.susceptibility_table.csv',header = TRUE, sep = ",")
n1=2500
n2=2500

M = matrix(, nrow = 100, ncol = 3)
for(i in 1:100){
	p = pow_func(S[i,])
	M[i,1] = S[i,1]
	M[i,2] = p[[1]]
	M[i,3] = p[[2]]
}

Mt = as_tibble(M)
write_csv(Mt,"../../examples/study_power.p_vs_s.csv")


# back in python
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns

f = open("../../examples/study_power.p_vs_s.csv")
header = f.readline()
X = [[float(x) for x in l.strip().split(',')] for l in f]
f.close()

X = np.array(X)
S = np.hstack([X[:,0],X[:,0]])
Ph = X[:,1]
Pl = X[:,2]
P = np.hstack([Ph,Pl])
mode = ['incidence_high']*X.shape[0] + ['incidence_low']*X.shape[0]
df = {'Susceptibility': S, 'Power': P, 'incidence model': mode}
df = pd.DataFrame(df)
ax = sns.lineplot(data=df, x='Susceptibility', y='Power', hue='incidence model')
_=plt.tight_layout()
plt.savefig('../../examples/study_power.p_vs_s.pdf')
plt.close()

