#import autograd.numpy as np 
import numpy as np
from scipy.special import logsumexp
import collections

def tups2dicts(Ll_tups):
    '''
    Take in a list(corresponding to one individual) of lists(corresponding to chromosomes) of tups and converts each nested element into a dict. Returns list with nesting structure intact
    '''
    outlist = []
    for chr_list_tups in Ll_tups:
        chr_list_dicts = [{'type':x[0],'length':x[1]} for x in chr_list_tups]
        outlist.append(chr_list_dicts)
    return(outlist)      


def computeLoglikelihood_binomial(D,P,alpha=14.3):
    '''
    Takes in data dictionary D ({'10':1018, '11':2227.16, '00':100.667}) and outputs loglikelihood (as well aslikelihood) under binomial model of diploid ancestry states. Uses alpha as a hyperparameter to account for admixture LD
    '''
#{{{    

    for key in ('10','11','00'): #Set key value to 0 if key not in dictionary
        if key not in D.keys():
            D[key]=0

    assert all([elem in ('10','11','00') for elem in D.keys()])
    p_A,p_B = P[0],P[1]
        
    hom1 = D['11']*(np.log(1-p_A) + np.log(1-p_B)) if D['11'] != 0 else 0
    hom2 = D['00']*(np.log(p_A) + np.log(p_B)) if D['00'] != 0 else 0

    het  = D['10']*_logsum(np.log(p_A) + np.log(1-p_B), np.log(p_B) + np.log(1-p_A)) if D['10'] !=0 else 0

    loglik = (1/alpha)*(hom1 + hom2  + het)  #Divide by 1/alpha to deal with non-independence - alpha depends on population

#}}}
    return(loglik)


"~~~~~~~~~~~~~~~~~~~~~~~"
def _logsum(logx, logy):
#{{{ 
    '''
    Evaluates log(x+y)
    '''
    if logx == -np.inf : 
        return(logy)
    if logy == -np.inf:
        return(logx)

    if logx > logy:
        return logx + np.log(1 + np.exp(logy-logx))
    else:
        return logy + np.log(1 + np.exp(logx-logy))

#}}}


def computeParams_markovWF(pars,e):
    """
    Applies the markovian wright-fisher based model to compute lambdaA's and lambdaB's taking into account hidden ancestry switches. Adds a small e to avoid 0 lambdas.
    """

    p_A,p_B,t_A,t_B = pars[0] , pars[1] , pars[2], pars[3]
    params = {'lambdaA':[(1-p_A)*t_A + e, p_A*t_A + e], 'lambdaB':[(1-p_B)*t_B + e, p_B*t_B + e]} 

    return(params)

def precompute_transitionMatrix(params,mode ='regular'):
    '''
    Takes in dict of parameter values ( e.g : {'lambdaA':[0.5,0,5], 'lambdaB':[0.3,0.4]})
    Returns precomputed parts of transition probabilites matrix N + E ( broken into non-exponential parts and exponential parts)
    '''
#{{{ 
    N = np.array([[0, np.log(params['lambdaA'][1]), np.log(params['lambdaB'][1]), -np.inf],
                 [np.log(params['lambdaA'][0]), 0,-np.inf ,np.log(params['lambdaB'][1])],
                 [np.log(params['lambdaB'][0]), -np.inf,0 ,np.log(params['lambdaA'][1])],
                 [-np.inf, np.log(params['lambdaB'][0]), np.log(params['lambdaA'][0]),0]
                ])


    E = np.array([[-(params['lambdaA'][0] + params['lambdaB'][0]), -(params['lambdaA'][1] + params['lambdaB'][0]), -(params['lambdaA'][0] + params['lambdaB'][1]),-np.inf],
                 [-(params['lambdaA'][0] + params['lambdaB'][0]),-(params['lambdaA'][1]+params['lambdaB'][0]),-np.inf,-(params['lambdaA'][1] + params['lambdaB'][1])],
                 [-(params['lambdaA'][0] + params['lambdaB'][0]), -np.inf,-(params['lambdaA'][0]+params['lambdaB'][1]),-(params['lambdaA'][1] + params['lambdaB'][1])],
                 [-np.inf,-(params['lambdaA'][1] + params['lambdaB'][0]),-(params['lambdaA'][0] + params['lambdaB'][1]) ,-(params['lambdaA'][1]+params['lambdaB'][1])]
                 ])

    if mode == 'cns': #Only non exponential part changes upon censoring. 
        N_cns = np.array([[-(np.log(params['lambdaA'][0]) + np.log(params['lambdaB'][0])), -np.log(params['lambdaB'][0]), -np.log(params['lambdaA'][0]), -np.inf],
                 [ -np.log(params['lambdaB'][0]), - (np.log(params['lambdaA'][1]) + np.log(params['lambdaB'][0])), -np.inf , -np.log(params['lambdaA'][1]) ],
                 [-np.log(params['lambdaA'][0]), -np.inf,-(np.log(params['lambdaA'][0]) + np.log(params['lambdaB'][1])) , -np.log(params['lambdaB'][1])],
                 [-np.inf, -np.log(params['lambdaA'][1]), -np.log(params['lambdaB'][1]),-(np.log(params['lambdaA'][1]) + np.log(params['lambdaB'][1]))]
                ])
        return(N,N_cns,E)
#}}}
    return(N,E)

    
def get_orderedSates(d):
    """
    Takes in a dict with a single key value pair. Uses hardcoded dictionaries to translate unordered states to ordered states.
    Returns ordered states corresponding to unordered state represented in d
    """
#{{{ 

    unordOrdDict_labs = {'10':[1,2], '11':[3,3], '00':[0,0]}
    unordOrdDict_states = {'00':[[0,0],[0,0]], '10':[[1,0],[0,1]] , '11':[[1,1],[1,1]]}

    unord_state = d['type']
    ord_labels = unordOrdDict_labs[unord_state]
    ord_states = unordOrdDict_states[unord_state]

#}}}
    return(ord_labels, ord_states)



def err_func(d,k,bg=None):
    '''
    d = tract dict with 'type' and 'length' , k = hyperparameter for exponential, bg = {'n': median number of tracts across chromosomes, 'type':background tract type}
    '''
#{{{ 

    x = d['length']/100
    
    if bg==None:
        loge = -k*x

    else:
        if (bg['n'] <= 5) and not (d['type'] == bg['type'] ):
            if 0 < x <= 0.05:
                loge = 0
            else:
                loge = -k*x
        else:
            loge = -k*x

#}}}
    return(loge)


def err_func_v0(d,k):
    x = d['length']/100
    return(np.exp(-k*x))

def get_emissions(d,d_prev,k,loge):
 #{{{   
    nonErrorLogProb=np.log(1-np.exp(loge))
    
    if k == 1 or k == 2:
        if d_prev['type'] == '00':
            emiss = [loge,nonErrorLogProb,-np.inf,-np.inf]
        elif d_prev['type'] == '11':
            emiss = [-np.inf,nonErrorLogProb,-np.inf,loge]
        else:
            emiss = [loge-np.log(2),nonErrorLogProb,-np.inf,loge-np.log(2)] #error probs are dsitributed across '00' and '11' if prev state is het

    elif k==0:
        if d_prev['type'] == '00':
            emiss = [nonErrorLogProb,loge-np.log(2),-np.inf,-np.inf]
        if d_prev['type'] == '0':
            emiss = [nonErrorLogProb,loge-np.log(2),-np.inf,-np.inf]

    else:
        sys.exit('invalid key')
#}}}
    return()

def computeLoglikelihood_cnsPM(D, pars, phi=69.314,err=False):

    '''
    Takes in dict of parameter values ( e.g : [0.3,0.5,6,5]}) and the data D for ALL chromsomes ({'type':'10', 'length':1018} ... ). Additionally takes in an error function and parameters phi for the error function ( see err_func). Also takes in bg which is a dictionary containing median number of tracts and background tract summary statistics from individual for use with error function. Computes loglikelihood of data under our pooled markov model with censoring ( cns). 
    '''
#{{{ 
    assert len(D) == 22
    
    #Compute background summary statistics across all chromosomes for use with error function
    types = [x['type'] for x in [x for y in D for x in y]]
    most_common = collections.Counter(types).most_common()[0][0]
    
    bg = {'n':np.median([len(x) for x in D]),'type':most_common}
    #Compute loglik across all chromosomes
    loglik=0

    print(err)
    if err:
        print("Running PAPI in error model mode")
        for d in D:
            loglik+=computeLoglikelihood_cnsPM_single_err(d,pars,phi,bg)
    else:
        print("Running PAPI in default (no error model) mode")
        for d in D:
            loglik+=computeLoglikelihood_cnsPM_single(d,pars)

#}}}
    return(loglik)


def computeLoglikelihood_cnsPM_single(D, pars):

    '''
    Takes in dict of parameter values ( e.g : [0.3,0.5,6,5]}) and the data D for a single chromosome ({'type':'10', 'length':1018} ... ). Additionally takes in an error function and parameters phi for the error function ( see err_func). Computes loglikelihood of data under our pooled markov model with censoring ( cns). 
    '''
#{{{ 

    #Get params (lambdas)
    params = computeParams_markovWF(pars,1e-6)
    assert all([ x != 0 for x in params.values()]) 
    #State space dict

    StSp = {0:[0,0],1:[1,0],2:[0,1],3:[1,1]}

    typDict = {'00':0,'10':1,'01':2,'11':3}
    typDict_phs = {'00':[[0,0],[0,0]],'10':[[1,0],[0,1]],'11':[[1,1],[1,1]]} #Used only when there is one censored tract

    
    #params = {'lambdaA':[(1-p_A)*t_A, p_A*t_A], 'lambdaB':[(1-p_B)*t_B, p_B*t_B]} 
    #Initialize P & Em matrices (2 * N )
    #Compute entries of P matrix (2 * N ) 
    if len(D) > 1:

        P = np.empty((4,len(D)))

        #Initialization
        for k in [0,1,2,3]:
            Em = {0:np.log([1, 0, 0 ,0]), 1:np.log([0,1,0,0]),2:np.log([0,1,0,0]),3:np.log([0,0,0,1])}  # Observed data is only '10'. 1 and 2 have identical emissions under this scheme
            P[k,0] = np.log(params['lambdaA'][StSp[k][0]]) - params['lambdaA'][StSp[k][0]]*(D[0]['length']/100)  + np.log(params['lambdaB'][StSp[k][1]]) - params['lambdaB'][StSp[k][1]]*(D[0]['length']/100) + Em[k][typDict[D[0]['type']]]

        N, N_cns, E = precompute_transitionMatrix(params,mode='cns')
        
        for col in range(1, len(D)): #Oth columns has already been computed 

            Em = {0:np.log([1, 0, 0 ,0]), 1:np.log([0,1,0,0]),2:np.log([0,1,0,0]),3:np.log([0,0,0,1])}

            if col == len(D)-1: #Censor last column
                T = N_cns + E*(D[col]['length']/100)
                #print('N_cns:\n {}'.format(N_cns))
                #print('E:\n {}'.format(E))
                #print('Transition Matrix :\n {}'.format(T))
            else:
                T = N + E*(D[col]['length']/100) #Divide by hundred to convert centimorgan to morgan

                #print('N:\n {}'.format(N))
                #print('E:\n {}'.format(E))
                #print('Transition Matrix :\n {}'.format(T))
            for k in [0,1,2,3]:
                P[k,col]= logsumexp([\
                        P[0,col-1] + T[0,k] + Em[k][typDict[D[col]['type']]],\
                        P[1,col-1] + T[1,k] + Em[k][typDict[D[col]['type']]],\
                        P[2,col-1] + T[2,k] + Em[k][typDict[D[col]['type']]],\
                        P[3,col-1] + T[3,k] + Em[k][typDict[D[col]['type']]],\
                        ])
            
    else: #Only one censored tract when len(D) == 1
        P = np.empty((2,len(D)))
        #print('One tract only , censoring')
        for k in [0,1]:
            phase = typDict_phs[D[0]['type']][k]
            P[k,0] = -( params['lambdaA'][phase[0]]*(D[0]['length']/100) + params['lambdaB'][phase[1]]*(D[0]['length']/100) ) 
   
    
#}}}
    return(logsumexp(P[:,-1]))



def computeLoglikelihood_cnsPM_single_err(D, pars, phi,bg):

    '''
    Takes in dict of parameter values ( e.g : [0.3,0.5,6,5]}) and the data D for a single chromosome ({'type':'10', 'length':1018} ... ). Additionally takes in an error function and parameters phi for the error function ( see err_func). Computes loglikelihood of data under our pooled markov model with censoring ( cns). 
    '''
#{{{ 

    #Get params (lambdas)
    params = computeParams_markovWF(pars,1e-6)
    assert all([ x != 0 for x in params.values()]) 
    #State space dict

    StSp = {0:[0,0],1:[1,0],2:[0,1],3:[1,1]}

    typDict = {'00':0,'10':1,'01':2,'11':3}
    typDict_phs = {'00':[[0,0],[0,0]],'10':[[1,0],[0,1]],'11':[[1,1],[1,1]]} #Used only when there is one censored tract

    
    #params = {'lambdaA':[(1-p_A)*t_A, p_A*t_A], 'lambdaB':[(1-p_B)*t_B, p_B*t_B]} 
    #Initialize P & Em matrices (2 * N )
    #Compute entries of P matrix (2 * N ) 
    if len(D) > 1:

        P = np.empty((4,len(D)))

        #Initialization
        loge = err_func(D[0], phi,bg)
        for k in [0,1,2,3]:

            Em = {0:[np.log(1-np.exp(loge)), loge - np.log(2), -np.inf ,loge-np.log(2)],\
                    1:[loge-np.log(2),np.log(1-np.exp(loge)),-np.inf,loge-np.log(2)],\
                    2:[loge-np.log(2),np.log(1-np.exp(loge)),-np.inf,loge-np.log(2)],\
                    3:[loge-np.log(2),loge-np.log(2),-np.inf,np.log(1-np.exp(loge))]\
                    }  # Observed data is only '10'. 1 and 2 have identical emissions under this scheme

            P[k,0] = np.log(params['lambdaA'][StSp[k][0]]) - params['lambdaA'][StSp[k][0]]*(D[0]['length']/100)  + np.log(params['lambdaB'][StSp[k][1]]) - params['lambdaB'][StSp[k][1]]*(D[0]['length']/100) + Em[k][typDict[D[0]['type']]]

        N, N_cns, E = precompute_transitionMatrix(params,mode='cns')
        
        for col in range(1, len(D)): #Oth columns has already been computed 
            loge = err_func(D[col], phi,bg)


            Em = {0:[np.log(1-np.exp(loge)), loge - np.log(2), -np.inf ,loge-np.log(2)],\
                    1:[loge-np.log(2),np.log(1-np.exp(loge)),-np.inf,loge-np.log(2)],\
                    2:[loge-np.log(2),np.log(1-np.exp(loge)),-np.inf,loge-np.log(2)],\
                    3:[loge-np.log(2),loge-np.log(2),-np.inf,np.log(1-np.exp(loge))]\
                    }  # Observed data is only '10'. 1 and 2 have identical emissions under this scheme

            if col == len(D)-1: #Censor last column
                T = N_cns + E*(D[col]['length']/100)
                #print('N_cns:\n {}'.format(N_cns))
                #print('E:\n {}'.format(E))
                #print('Transition Matrix :\n {}'.format(T))
            else:
                T = N + E*(D[col]['length']/100) #Divide by hundred to convert centimorgan to morgan

                #print('N:\n {}'.format(N))
                #print('E:\n {}'.format(E))
                #print('Transition Matrix :\n {}'.format(T))
            for k in [0,1,2,3]:
                P[k,col]= logsumexp([\
                        P[0,col-1] + T[0,k] + Em[k][typDict[D[col]['type']]],\
                        P[1,col-1] + T[1,k] + Em[k][typDict[D[col]['type']]],\
                        P[2,col-1] + T[2,k] + Em[k][typDict[D[col]['type']]],\
                        P[3,col-1] + T[3,k] + Em[k][typDict[D[col]['type']]],\
                        ])
            
    else: #Only one censored tract when len(D) == 1
        P = np.empty((2,len(D)))
        #print('One tract only , censoring')
        for k in [0,1]:
            phase = typDict_phs[D[0]['type']][k]
            P[k,0] = -( params['lambdaA'][phase[0]]*(D[0]['length']/100) + params['lambdaB'][phase[1]]*(D[0]['length']/100) ) 
   
    
#}}}
    return(logsumexp(P[:,-1]))


#def computeLoglikelihood_cnsPM_flaterr(D, pars,e):
#    '''
#    Takes in dict of parameter values ( e.g : [0.3,0.5,6,5]}) and the data D for a single chromosome ({'type':'10', 'length':1018} ... ). Computes loglikelihood of data under our pooled markov model with censoring ( cns)
#    '''
##{{{ 
#
#    #Get params (lambdas)
#    params = computeParams_markovWF(pars,1e-6)
#    print(params)
#    assert all([ x != 0 for x in params.values()]) 
#    #State space dict
#
#    StSp = {0:[0,0],1:[1,0],2:[0,1],3:[1,1]}
#
#    typDict = {'00':0,'10':1,'01':2,'11':3}
#    typDict_phs = {'00':[[0,0],[0,0]],'10':[[1,0],[0,1]],'11':[[1,1],[1,1]]}
#
#    Em = {0:np.log([1-e, e/2, 0 ,e/2]), 1:np.log([e/2,1-e,0,e/2]),2:np.log([e/2,1-e,0,e/2]),3:np.log([e/2,e/2,0,1-e])}  # Observed data is only '10'. 1 and 2 have identical emissions under this scheme
#
#    #params = {'lambdaA':[(1-p_A)*t_A, p_A*t_A], 'lambdaB':[(1-p_B)*t_B, p_B*t_B]} 
#    #Initialize P & Em matrices (2 * N )
#    #Compute entries of P matrix (2 * N ) 
#    if len(D) > 1:
#
#        P = np.empty((4,len(D)))
#
#        #Initialization
#        for k in [0,1,2,3]:
#
#            P[k,0] = np.log(params['lambdaA'][StSp[k][0]]) - params['lambdaA'][StSp[k][0]]*(D[0]['length']/100)  + np.log(params['lambdaB'][StSp[k][1]]) - params['lambdaB'][StSp[k][1]]*(D[0]['length']/100) + Em[k][typDict[D[0]['type']]]
#
#        N, N_cns, E = precompute_transitionMatrix(params,mode='cns')
#
#        for col in range(1, len(D)): #Oth columns has already been computed 
#
#            if col == len(D)-1: #Censor last column
#                print('censoring')
#                T = N_cns + E*(D[col]['length']/100)
#                #print('N_cns:\n {}'.format(N_cns))
#                #print('E:\n {}'.format(E))
#                #print('Transition Matrix :\n {}'.format(T))
#            else:
#                T = N + E*(D[col]['length']/100) #Divide by hundred to convert centimorgan to morgan
#
#                #print('N:\n {}'.format(N))
#                #print('E:\n {}'.format(E))
#                #print('Transition Matrix :\n {}'.format(T))
#            for k in [0,1,2,3]:
#                P[k,col]= logsumexp([\
#                        P[0,col-1] + T[0,k] + Em[k][typDict[D[col]['type']]],\
#                        P[1,col-1] + T[1,k] + Em[k][typDict[D[col]['type']]],\
#                        P[2,col-1] + T[2,k] + Em[k][typDict[D[col]['type']]],\
#                        P[3,col-1] + T[3,k] + Em[k][typDict[D[col]['type']]],\
#                        ])
#
#    else: #Only one censored tract when len(D) == 1
#        P = np.empty((2,len(D)))
#        print('One tract only , censoring')
#        for k in [0,1]:
#            phase = typDict_phs[D[0]['type']][k]
#            P[k,0] = -( params['lambdaA'][phase[0]]*(D[0]['length']/100) + params['lambdaB'][phase[1]]*(D[0]['length']/100) ) 
#    
#
#    #print('Final P:{}'.format(P))
##}}}
#    return(logsumexp(P[:,-1]))



