# run-analytic-models
Script that can be run from the CLI to generate radio galaxy dynamical models from https://github.com/mhardcastle/analytic. This also enables the user to provide a set of observational constraints which can be fitted to the models to determine a set of radio galaxy models that describe a particular source.

#run instructions
-Make a directory to store all data
-In that directory clone the analytic models repo from https://github.com/mhardcastle/analytic
-In the same directory clone this repo
-In the same directory make a folder for your runs
-In that directory run:

python ../<rundynmodelsscriptdir>/run_dynmodels.py

This runs dynmaical models with a set of default parameters

If you want to run with your own parameters, then specify those options on the CLI:

python ../<rundynmodelsscriptdir>/run_dynmodels.py -Qjet "[1e37,1e39]" -m500 "[1e13,1e14,2e14]" -zeta "[0.5,1.0]"

This will run the analytic model for each of the specified jet powers, m500 and zeta (fraction of energy in B-field to radiating particles, calculate their synchrotron emission and save the outputs to a file 'output.pickle'


