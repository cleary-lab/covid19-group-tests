# Generating power curves for IgG study
All scripts to be run in the current directory

## Simulate an epidemic process
Use the scripts in study_power.epidemic_process.R to generate infection incidence curves over time. For each curve, choose an R0 and a prevalence of already infected individuals ("recovered"). It could also make sense to vary the initial incidence of the infected population on day 0. Change the length of the epidemic simulation if desired (default=365 days).

This will generate an SEIR plot for each model, and save incidence values in csvs to be used in the next step.

## Calculate the expected proportion of NP+ results
Run from the command line (requires numpy), e.g.:
```bash
python study_power.susceptibility_probs.py --incidence-high ../../examples/study_power.incidence.high.csv --incidence-low ../../examples/study_power.incidence.low.csv --sampling-times 0,30,60,90,120,150,180 --outfile ../../examples/study_power.susceptibility_table.csv
```

## Calculate power in the fisher exact test
Use the scripts in study_power.fisher_power.R to select a significance level (default=0.05) and the initial size of IgG- / IgG+ populations, and calculate power. If you have basic R skills, you can directly plot the results from here. Otherwise, save the results and use the python script below to plot.

## Plot the results
Run from the command line (requires matplotlib and seaborn):
```bash
python study_power.plot.py --infile ../../examples/study_power.p_vs_s.csv --infile ../../examples/study_power.p_vs_s.pdf
```
