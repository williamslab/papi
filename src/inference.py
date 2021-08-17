import pymc3 as pm
import models as md
import arviz as az
import theano.tensor as tt
import pandas as pd
import numpy as np
import pickle
import json
import sys
import scipy
from scipy import optimize
# define a theano Op for our likelihood function

class LogLike(tt.Op):

    """
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike, data):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        Parameters
        ----------
        loglike:
            The log-likelihood (or whatever) function we've defined
        data:
            The "observed" data that our log-likelihood function takes in
        """

        # add inputs as class attributes
        self.likelihood = loglike
        self.data = data

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        theta, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta,self.data)

        outputs[0][0] = np.array(logl) # output the log-likelihood

def filter_tracts(llt):
    '''
    Input list of lists of tups (llt) , Process so that len 0 tracts are dropped, and convert '01''s to '10's
    '''
    llt_proc=[]
    for lt in llt:
        lt_proc_zeros = [x for x in lt if not x[1]==0]
        lt_proc_zeros_het =[('10',x[1]) if x[0]=='01' else x for x in lt_proc_zeros]
        assert all([x[1]>0 for x in lt_proc_zeros_het]),"negative lengths in tract input"

        llt_proc.append(lt_proc_zeros_het)

    
    return(llt_proc)

def lik_func(params,D):
    
    D_dicts = md.tups2dicts(D)
    D_flat = [{x['type']:x['length']} for y in D_dicts for x in y]
    counter = collections.Counter() 
    for D_i in D_flat:  
        counter.update(D_i)     
    D_flat = dict(counter)
    l = md.computeLoglikelihood_binomial(D_flat,params[:2]) + md.computeLoglikelihood_cnsPM(D_dicts,params)
    #l = md.computeLoglikelihood_cnsPM(D_dicts,params)
    return(l)


def estimate_MAP(d_tracts):
    '''
    Function that infers map and CIs using the scipy.optimize module 
    '''
    
    bnds = ((0.001, 0.999),(0.001,0.999),(0,25), (0, 25))
    res = scipy.optimize.minimize(
        fun=lambda params, D: -lik_func(params,d_tracts),
        x0=np.array([0.7,0.8,4,7]),
        args=(d_tracts,),
        method='L-BFGS-B',
        bounds=bnds
    )
    
    return(np.round(res.x,2))


if __name__ == "__main__":
#{{{ 
    import ast 
    import collections 
    import numpy as np 
    import argparse

    # initiate the parser
    parser = argparse.ArgumentParser()

    # add long and short argument
    parser.add_argument("--inputfile", "-i",required=True, help="set input file",type=str)
    parser.add_argument("--ind", "-ind", type=int,required=True, help="set individual on which to run hmm on")
    parser.add_argument("--tracefile", "-t", nargs='?', const=None)
    parser.add_argument("--outfile", "-o",type=str, required=True)
    parser.add_argument("--mode", "-m", nargs='?', type=str, const='pymc')

    # read arguments from the command line
    args = parser.parse_args()

    #Read in data and header 
    data_all = pd.read_table(args.inputfile)
    data_ind = data_all.iloc[args.ind-1] 
    d_tracts = filter_tracts(ast.literal_eval(data_ind.tracts))


    if args.mode==None or args.mode == 'pymc':

        print("Running inference in pymc mode")
        basic_model = pm.Model()
        # create our Op
        logl = LogLike(lik_func, d_tracts)
        
        # use PyMC3 to sampler from log-likelihood
        with basic_model:
            # Priors for unknown model parameters
            p1=  pm.Uniform("p1",0,1)
            p2 = pm.Uniform("p2", 0, 1)
            t1 = pm.Uniform("t1", 0, 25)
            t2 = pm.Uniform("t2", 0, 25)
        
            # convert m and c to a tensor vector
            theta = tt.as_tensor_variable([p1,p2,t1,t2])
        
            # use a DensityDist (use a lamdba function to "call" the Op)
            pm.DensityDist('likelihood', lambda v: logl(v), observed={'v': theta})

            mean_q = pm.find_MAP(model=basic_model)
            #std_q = ((1/pm.find_hessian(mean_q))**0.5)[0]

            if args.tracefile:
                #full mcmc
                print("Running Full MCMC")
                trace = pm.sample(1500)
                df = az.summary(trace,credible_interval=0.9)
                df.to_csv(args.outfile[:-4] + '.summ',sep='\t')
                pm.save_trace(trace, args.tracefile,overwrite=True) 
                
            
            with open(args.outfile, 'w') as f:
                print(mean_q, file=f)



    if args.mode=='scipy-optimize':
        print("Running inference in scipy-optimize mode")
        MAPestimate_scipy = estimate_MAP(d_tracts)
        np.savetxt(args.outfile, MAPestimate_scipy, newline=" ")

#}}}

