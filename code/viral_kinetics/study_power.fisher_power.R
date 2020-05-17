library(statmod)

sig_level <- 0.05

pow_func = function(pars){
p1 = power.fisher.test(pars$IgG_pos.NP_pos.high, pars$IgG_neg.NP_pos.high, n1, n2, alpha=sig_level, nsim=200, alternative="two.sided")
p2 = power.fisher.test(pars$IgG_pos.NP_pos.low, pars$IgG_neg.NP_pos.low, n1, n2, alpha=sig_level, nsim=200, alternative="two.sided")
return(list(p1,p2))
}
S = read.table('../../examples/study_power.susceptibility_table.csv',header = TRUE, sep = ",")

# Initial population of IgG+
n1=2500
# Initial population of IgG-
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