# back in python
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--infile', help='Power vs susceptibility table, e.g. ../../examples/study_power.p_vs_s.csv')
	parser.add_argument('--outfile', help='Filename of saved plot, e.g. ../../examples/study_power.p_vs_s.pdf')
	args,_ = parser.parse_known_args()
	f = open(args.infile)
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
	plt.savefig(args.outfile)
	plt.close()