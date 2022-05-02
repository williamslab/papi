import numpy as np
import random 
import itertools
import sys 
import re 

def closest(lst, K): 
     lst = np.asarray(lst) 
     idx = (np.abs(lst - K)).argmin() 
     return lst[idx]

def any_discordantPairs(s):
    '''
    Checks if there are any discordant pairs ( non overlapping pairs) in the string
    E.g :
    '10112233' returns True while '11223344' returns False
    '''
    
    pairs = [s[index]+s[index+1] for index,elem in enumerate(s) if not index%2]
    answer = any([x != len(x) * x[0] for x in pairs])
    
    return(answer)

    
def generateAdmixStrings_branch1(t):
    '''
    Generates all possible admixture strings for founders on one branch given t.
    Prunes out cases without at least one '01' pair
    '''
    #{{{
    t = t+1
    base_case = ['1','0'] #Admixture strings for a pair of founders. This is also the resutl for the trivial case of t=1
    if t == 1:
        return(base_case)
    admixture_strings = [base_case]

    
    for iterations in range(t-1):
        current_branch_strings = admixture_strings[iterations] 
        next_branch_strings = []
    
        for i in range(len(current_branch_strings)): #Loop over the strings and stich together all possible unique combinations
            for j in range(i,len(current_branch_strings)): 
                next_branch_strings.append(current_branch_strings[i] + current_branch_strings[j])
        
        if iterations == 0: #remove '00' from t=2 case
            next_branch_strings.remove('00')
        admixture_strings.append(next_branch_strings)

    #}}}
    return(admixture_strings[-1])


def generateAdmixStrings1(t1,t2):
    '''
    Function that generates admixture strings for both branches with admixture occuring t1 and t2 generations ago respectively 
    '''
    #{{{ 
    returnStrings = []
    if t1 == t2:
        branch_strings = generateAdmixStrings_branch(t1)
        
        for i in range(len(branch_strings)):
            for j in range(i,len(branch_strings)):
                returnStrings.append(branch_strings[i] + branch_strings[j])
        
        if t1 == 0:
            returnStrings.remove('00')
    else:
        t_max = max(t1,t2)
        t_min = min(t1,t2) 
        diff = t_max - t_min

        branch_strings_min = [''.join([char*(2**diff) for char in s]) for s in generateAdmixStrings_branch(t_min)] #Duplicates every string in a list of strings ( e.g '101' become '110011'). Duplication extent is based on diff between t's
        branch_strings_max = generateAdmixStrings_branch(t_max)
        
        for elem_max in branch_strings_max:
            for elem_min in branch_strings_min:
                returnStrings.append(elem_max + elem_min)
    #}}}
    return(returnStrings)

def label_string(s):
    '''
    Outputs labelled string according to map_dict ( dict should map pairs of characters onto single character )
    '''
    map_dict = {'10':'A', '11':'1', '00':'0'}
    pairs = [s[index]+s[index+1] for index,elem in enumerate(s) if not index%2]
    mapped_string = ''.join([map_dict[pair] for pair in pairs])

    return(mapped_string)

def label_string_v2(s, t1, t2):
    '''
    Outputs labelled string according to difference between t's to provide the shortest representation possible
    '''
#{{{ 
    
    #Label first half corresponding to larger t 
    s_first_half = s[0:len(s)//2]
    s_first_half_shortened = label_string(s_first_half)

    
    #Lable second half corresponding to shorter t ( successive labelling depends on the difference between the two ts)
    s_second_half = s[len(s)//2:len(s)]
    diff = max(t1,t2) - min(t1,t2)
    
    s_second_half_shortened = s_second_half
    for iteration in range(diff):
        map_dict = {'11':'1', '00':'0'}
        pairs = [s_second_half_shortened[index]+s_second_half_shortened[index+1] for index,elem in enumerate(s_second_half_shortened) if not index%2]
        s_second_half_shortened = ''.join([map_dict[pair] for pair in pairs])

    mapped_string = s_first_half_shortened + s_second_half_shortened
#}}}
    return(mapped_string)

def bucketAdmixstrings(admixstrings):
    '''
    Puts admixstrings into a dictionary with key representing p1/p2 - e.g '20/50' - and entry is a list of admixstrings 
    '''
    #{{{ 
    plist = np.array([0.1,0.2,0.3,0.4,0.5])
    returnDict = dict()

    for s in admixstrings:
        s_l = s[:len(s)//2]
        s_r = s[len(s)//2:]

        pl = s_l.count('0')/len(s_l)
        pr = s_r.count('0')/len(s_r)
        
        if (pl != 0) and (pl != 1) :
            pl_match = plist[np.abs(plist-pl).argmin()]
        else:
            pl_match = pl
    
        if (pr != 0) and (pr != 1) :
            pr_match = plist[np.abs(plist-pr).argmin()]
        else:
            pr_match = pr

        key = str(pl_match*100)[:-2] + '_' +  str(pr_match*100)[:-2]
        if key in returnDict:
            returnDict[key].append(s)
        else:
            returnDict[key] = [s] 
#}}}
    return(returnDict)

############# String generators(linear growth) ###############

def translateString(s):
    trans_dict = {'A':'10','1':'11','0':'00'}
    translated_string=''.join([trans_dict[x] for x in s])
    return(translated_string)

def generateAdmixStrings_branch(t,p,n=2):
    '''
    generate shuffled admixstrings for one branch given one value of p. default number of shuffles for each string is 2. No subadmixture
    '''
#{{{ 
    if t==0:
        return(['0','1'])

    else:
        num_founder_branches = 2**(t-1)
        
        if p<=0.5:
            base_string = '1'*num_founder_branches
        else:
            p = 1-p
            base_string = '0'*num_founder_branches

        prop_admixed_branches = p*2
        allowed_prop_admixed_branches = [x/num_founder_branches for x in range(num_founder_branches+1)]
	    
        num_admixed_branches = int(closest(allowed_prop_admixed_branches,prop_admixed_branches)*num_founder_branches)
        
        rem_shufs=n 

        strings = []
        while rem_shufs > 0:
            replace_inds = random.sample(range(0, len(base_string)), num_admixed_branches)
            rem_shufs-=1
            string=base_string

            for i in replace_inds: #Replace sampled founder branches with 'A' for admixed
                string=string[:i] +'A' + string[i+1:]

            
            final_string=translateString(string) #Translate branches string to admixstring
            strings.append(final_string)
#}}}
    return(list(set(strings)))

def generateAdmixStrings_branch_subadmix(t,p,nl=2,opts='default'):
    '''
    generate admixstrings for on branch with subadmixture given t and p. argument nl is number of different permutations to make for a given number of admixed founder (sub)branches on the left branch in the subadmixture. Right sub-branch only gets one permutation per number of admixed founder (sub)branches.
    '''
#{{{ 
    assert t>=3

    tl=t-1 
    strings_l = generateAdmixStrings_branch(tl,p,nl)
    strings_r_nested=[]

    if opts == 'all':
        subadmix_trs = list(reversed(range(1,t)))
    elif opts == 'default':
        subadmix_trs = list(set([t//2,1]))
    else:
        print('invalid option : {}'.format(opts))
        sys.exit()


    for tr in subadmix_trs:
        strings_r = generateAdmixStrings_branch(tr,p,1)
        diff = tl-tr
        strings_r_adjusted = [''.join([char*(2**diff) for char in s]) for s in strings_r] 
        strings_r_nested.append(strings_r_adjusted)


    strings_r= [x for y in strings_r_nested for x in y]
    comb_list = list(itertools.product(strings_l, strings_r))
    return_list = list(set([x[0] + x[1] for x in comb_list]))
#}}}
    return(return_list)

def generateAdmixStrings(t1,t2,p1,p2,n=2,options=['n','n']):
    '''
    Function that generates admixture strings for both branches with admixture occuring t1 and t2 generations ago respectively. Option list contains 'n' for normal and 's' for subadmix. There's one option for each branch
    '''
    #{{{ 

    if options[0] == 'n':
        strings_l = generateAdmixStrings_branch(t1,p1,n) 
    elif options[0] == 's':
        strings_l = generateAdmixStrings_branch_subadmix(t1,p1,n)
    else:
        print('invalid option for t1')
        return()
    
    if options[0] == 'n':
        strings_r = generateAdmixStrings_branch(t2,p2,n) 
    elif options[0] == 's':
        strings_r = generateAdmixStrings_branch_subadmix(t2,p2,n)
    else:
        print('invalid option for t2')
        return()
    
    diff = t1-t2
    strings_r_adjusted = [''.join([char*(2**diff) for char in s]) for s in strings_r] 
    comb_list = list(itertools.product(strings_l, strings_r_adjusted))
    return_list = list(set([x[0] + x[1] for x in comb_list]))

    #}}}
    return(return_list)

def lambdaRecursive(admixStr):
    '''
    Recursive function that computes lambda for the transmitted focal haplotype given the founder admixture string on one pedigree branch
    '''
   #{{{ 
    #Base case are the transmitted haplotypes from the founder generation
    ps = [0 if s=='1' else 1 for s in admixStr]
    lambdas=[0]*(len(admixStr))


    het = lambda ptup: ptup[0]*(1-ptup[1]) + ptup[1]*(1-ptup[0])
    while len(ps) > 1:

        p_pairs = [(ps[index],ps[index+1]) for index,elem in enumerate(ps) if not index%2]
        lambda_pairs = [(lambdas[index],lambdas[index+1]) for index,elem in enumerate(lambdas) if not index%2]
        assert len(p_pairs) == len(lambda_pairs) == len(ps)//2


        #Get ps and lambdas for next gen using pairs
        ps = [(tup[0] + tup[1])/2 for tup in p_pairs]
        lambdas = [(l[0]+l[1])/2 + het(p_pairs[index]) for index,l in enumerate(lambda_pairs)]
        #print(ps,lambdas,[het(x) for x in p_pairs])
    #}}}
    return(lambdas)

def rate_plus_het(rate,het):


    if not rate=='99':
        result=rate+het
    else:
        if het==0:
            result=rate
        else:
            result=het

    return(result)

def lambdaRecursive_pops(admixStr):
    '''
    Recursive function that computes lambda for the transmitted focal haplotype given the founder admixture string on one pedigree branch
    '''
   #{{{ 
    #Base case are the transmitted haplotypes from the founder generation
    ps = [0 if s=='1' else 1 for s in admixStr]
    lambdas=[0 if s=='0' else '99' for s in admixStr] #tracks lambda_0 ( rate for label '0')


    het = lambda ptup: ptup[0]*(1-ptup[1]) + ptup[1]*(1-ptup[0])
    rate_add = lambda ltup: ltup[0] if ltup[1]=='99' else ltup[1] if ltup[0]=='99' else '99' if ltup[1]==ltup[0]=='99' else (ltup[0]+ltup[1])/2

    while len(ps) > 1:

        p_pairs = [(ps[index],ps[index+1]) for index,elem in enumerate(ps) if not index%2]
        lambda_pairs = [(lambdas[index],lambdas[index+1]) for index,elem in enumerate(lambdas) if not index%2]
        assert len(p_pairs) == len(lambda_pairs) == len(ps)//2

        #Get ps and lambdas for next gen using pairs
        ps = [(tup[0] + tup[1])/2 for tup in p_pairs]
        lambdas = [rate_plus_het(rate_add(l), het(p_pairs[index])) for index,l in enumerate(lambda_pairs)]
        #print(ps,lambdas,[het(x) for x in p_pairs])
    #}}}
    return(lambdas)

def displace_string(s,posList):
    '''
    Displace the second half of string s into the first half in units of '{11..}{00..}' given positions in the first half string.
    '''
    sL,sR = s[:len(s)//2],s[len(s)//2:]
    nL,nR=len(max(re.compile("1*").findall(sL)))*2,len(max(re.compile("1*").findall(sR)))*2
    split_sL,split_sR=[(sL[i:i+nL]) for i in range(0, len(sL), nL)],[(sR[i:i+nR]) for i in range(0, len(sR), nR)]
    
    for i in posList:
        split_sL.insert(i,split_sR.pop(0))
    
    return(("".join(split_sL + split_sR)))

