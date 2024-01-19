import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
import xarray as xr
import os
import sys
import time
import pickle
import scipy as sp
from scipy import stats
import math
import ast
label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size


def IntersecOfSets(arr1, arr2, arr3): 
    # Converting the arrays into sets 
    s1 = set(arr1) 
    s2 = set(arr2)
    if arr3 is None:
        s3 = []
    else: 
        s3 = set(arr3) 
      
    # Calculates intersection of  
    # sets on s1 and s2 
    set1 = s1.intersection(s2)         #[80, 20, 100] 
    if s3 == []:
        result_set = set1
    else:  
    # Calculates intersection of sets 
    # on set1 and s3 
        result_set = set1.intersection(s3) 
      
    # Converts resulting set to list 
    final_list = list(result_set) 
    return final_list

class log_uniform():        
    def __init__(self, a=-1, b=0, base=10):
        self.loc = a
        self.scale = b - a
        self.base = base

    def rvs(self, size=None, random_state=None):
        uniform = sp.stats.uniform(loc=self.loc, scale=self.scale)
        if size is None:
            return np.power(self.base, uniform.rvs(random_state=random_state))
        else:
            return np.power(self.base, uniform.rvs(size=size, random_state=random_state))

def arg_as_list(s):
    v = ast.literal_eval(s)
    if type(v) is not list:
        raise argparse.ArgumentTypeError("Argument \"%s\" is not a list" % (s))
    return v

def vals_as_list(x):
    #want to convert dict values into lists to allow recursive selection to find combinations
    for key, value in x.items():
        if not isinstance(value,list):
            x[key]=[value]
    return x


def float_or_floatlist(arg):
    if arg is None:
        return arg
    else:
        return np.float(arg)

def run_dynmodel(dict_):
    #import solver scripts to current working dir
    #sys.path.append(args['solverpath'][0])
    #os.system('cp ' + args['solverpath'][0] + '/*.py .')
    from solver import Evolve_RG
    from synch_constants import *
    #from solver import Evolve_RG
    #from constants import *
    print 'Running Evolve_RG with the following:\n'
    print dict_
    model=Evolve_RG(dict_['environment'],kT=dict_['kt'],p0=dict_['p0'],rc=dict_['rc'],beta=dict_['betain'],
        do_adiabatic=dict_['adiabaticlosses'],q=dict_['q'],z=dict_['z'],verbose=False,zeta=dict_['zeta'],gmin=dict_['gmin'],
        gmax=dict_['gmax'],xi=dict_['lobe_shock_ratio'],M500=dict_['m500'])

    #set the time steps
    tv = np.concatenate((np.logspace(-5,np.log10(dict_['tmin']),50)*Myr,np.linspace(1.01,dict_['tmax'],100)*Myr))
    print 'Now solving for ',dict_['Qjet']
    model.solve(dict_['Qjet'],tv)


    model.save(dict_['save'])

def main():
    
    import argparse
    import ast
    from sklearn.model_selection import ParameterGrid

    parser = argparse.ArgumentParser(description='Run semi-analytic models from Hardcastle+18 and fit to obsevation data')
    
    timeparser = parser.add_argument_group("-------------------------Time Settings-------------------------")
    timeparser.add_argument('-tmin', help="Start time (Myr)", default=0.1, type=int)
    timeparser.add_argument('-tmax', help="End time (Myr)", default=100, type=int)
    timeparser.add_argument('-tsteps', help="Time increment list (Myr)", default=100, type=list)

    envparser = parser.add_argument_group("-------------------------Environment Settings-------------------------")
    envparser.add_argument('-environment', help="Type of environment ('beta' or 'universal')", default='universal', type=str)
    envparser.add_argument('-m500', help="M500 mass in solar masses", default=[1e13,1e14], type=arg_as_list)
    envparser.add_argument('-betain', help="Beta index", default=0.7, type=arg_as_list)
    envparser.add_argument('-kt', help="Temperative kT in keV", default=5., type=float_or_floatlist)
    envparser.add_argument('-p0', help="Central pressure in Pa", default=1e-11, type=float_or_floatlist)
    envparser.add_argument('-rc', help="Core radius in kpc", default=100., type=float_or_floatlist)

    rgparser = parser.add_argument_group("-------------------------Radio Galaxy Settings-------------------------")
    rgparser.add_argument('-z', help="Redshift", default=0,type=float_or_floatlist)
    rgparser.add_argument('-q', help="Injection index", default=2.0,type=float_or_floatlist)
    rgparser.add_argument('-adiabaticlosses', help="Do adiabatic losses?", default=True,action='store_true')
    rgparser.add_argument('-zeta', help="Fraction of energy in radiating particles and B-field", default=1.,type=arg_as_list)
    rgparser.add_argument('-kappa', help="Fraction of non-radiating to radiating particles", default=0.,type=float_or_floatlist)
    rgparser.add_argument('-gmin', help="Minimum Lorentz factor of lobe electrons", default=10.,type=float_or_floatlist)
    rgparser.add_argument('-gmax', help="Maximum Lorentz factor of lobe electrons", default=1e6,type=float_or_floatlist)
    rgparser.add_argument('-lobe-shock-ratio', help="Ratio of energy in lobe to shock", default=0.5,type=float_or_floatlist)

    solveparser = parser.add_argument_group("-------------------------Solver Settings-------------------------")
    solveparser.add_argument('-Qjet', help='Jet power (W)',default=[10e38,10e39], type=arg_as_list)
    solveparser.add_argument('-findsynch', help='Calculate synchrotron emission',default=150e6 ,type=float_or_floatlist)
    solveparser.add_argument('-findcorrection', help='Calculate corrections for synchrotron emission',default=150e6 ,type=float_or_floatlist)

    saveparser = parser.add_argument_group("-------------------------Save Settings-------------------------")
    saveparser.add_argument('-save', help='Store solutions',default=None, type=str)

    parser.add_argument('--solverpath', help='Path to location of analytic directory', default='/beegfs/general/mahatmav/dynamicalmodelling/analytic/', type=str)


    options = parser.parse_args()
    args = vars(options)

    #determine combinations of model parameters
    listargs=vals_as_list(args)
    param_grid=ParameterGrid(listargs)

    #import solver scripts to current working dir
    sys.path.append(args['solverpath'][0])
    #os.system('cp ' + args['solverpath'][0] + '/*.py .')
    from solver import Evolve_RG
    #from constants import *

    #Start running the models
    for modelargs in param_grid:
        run_dynmodel(modelargs)


if __name__ == "__main__":
   main()

