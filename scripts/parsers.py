#Import required modules
import re   
import pandas as pd 
import numpy as np 

########### CORE FUNCTIONS ##############################

def haplo_parse(input_string):
    '''
    Args:
    input_string : Takes input string in bp format (output from ped-sim) 
    collapse : Collapses identical adjacent blocks by default
    drop_ends : drops final segment to avoid chromosome end bias. False by default
    Returns:
    haplo_dicts : list of dictionaries containing chromosome number, haplotype,
                  block length and [start, stop] positions as keys 

    '''
 #{{{
    
    input_string = input_string.rstrip() + ' ' #Remove newline and add space for better regex processing 
    parsed_list = re.findall(r'\d+\|.+?(?=\d{1,2}\||$)', input_string, re.MULTILINE) #Killer regex here. Captures information chromosome by chromosome 
    haplo_lists_of_dicts = [None]*len(parsed_list)

    
    
    for i in range(len(parsed_list)): #Loop over the chromosomes

        #regex parse each list entry to pull out all the information we need 
        chrom = re.search(r'^[^|]*',parsed_list[i]).group(0)
        chrom_start = float(re.search(r'\|(.*?)\s', parsed_list[i]).group(1))
        haplotypes = re.findall(r'\s(.*?):', parsed_list[i]) 
        chrom_end_list = re.findall(r':(.*?)\s', parsed_list[i])
        
        range_list = [[float(chrom_start),float(chrom_end_list[i])]  if i==0 else [float(chrom_end_list[i-1]),float(chrom_end_list[i])] for i in range(len(chrom_end_list))]
        end_status = [True if j == len(haplotypes) - 1 else False for j in range(len(haplotypes))]
        chrom_list = [{'haplotypes':haplotypes[j] ,'startstop':range_list[j], 'chromosome':chrom,\
                'length':range_list[j][1] - range_list[j][0], 'end_status':end_status[j] } \
                for j in range(len(haplotypes)) ] 

        haplo_lists_of_dicts[i] = chrom_list # Assing chromosome i's list of haplotype dictionaries to a list 

    
    
    haplo_dicts = [item for sublist in haplo_lists_of_dicts for item in sublist] #Collapse list 

    #If option collapse is provided ( provided by default) we collapse contiguous blocks and return list 
    haplo_dicts = collapse_blocks(haplo_dicts)

    #}}}
    return(haplo_dicts) 



def ancestryParse(input_strings, admixed_branches):
    '''
    Args:
    input_strings : Takes input strings - one for each chromosome - in bp format (output from ped-sim) 
    admixed_branches : list of admixture types for each branch ( list length should be equal to number of founders)

    Returns:
    output_strings : list of strings (one element for each chromosome) in bp format with haplotype id replaced by population/ancestry ids 
    '''
#{{{
    assert len(input_strings) == 2 

    #Get founder haplotypes using modular arithmetic ( based on number of meises simulated and current simulation number ) 
    sim_no = int(re.search('(\d+)_',input_strings[0]).group(1))
    num_founders = 2**(int(input_strings[0][0])+1) #correction of +1 because we use number of observed crossovers
    num_haplotypes = 2*num_founders
    
    assert len(admixed_branches) == num_founders # Check that both input strings correspond to admixed_branches input by checking how many meoises simulated ( 1st character of bp file will give us this information ) 


    first_id = (sim_no - 1)*num_haplotypes
    last_id = sim_no*num_haplotypes - 1 
    
    founder_haplotypes = [str(x) for x in range(first_id, last_id+1)]

    #Replace input string haplotypes with population identifiers based on admixed_branches
    ##Create dictionary mapping haplotypes to admixture types 
    admixed_branches_expanded = [admixed_branches[i//2] for i in range(len(admixed_branches)*2)] #Duplicate input for mapping to haplotypes in bp file 
    
    ##Replace input strings with corresponding admixture type
    output_strings = [None]*2
    for j in range(len(input_strings)): 
        output_strings[j] = input_strings[j]
        for i in range(len(founder_haplotypes)): 
            output_strings[j] = output_strings[j].replace(' ' + str(founder_haplotypes[i]) + ':', ' ' + admixed_branches_expanded[i] + ':') #Extra space is to ensure exact matches ( e.g, so that 16: doesn't match 6: )
            #}}}
    return(output_strings)


def admixSimu_parse(lines):
    '''
    Creates admixSimu parsed bp format lines from input bp lines by removing the first sequence of digits(ususally 0.0) after the '|'.
    '''

    newlines = []
    for line in lines:
        newline = re.sub(r'(\|)\d+', r'\g<1>', line)
        newlines.append(newline)

    return(newlines)

def computeAncprops(haplodicts):
    '''
    Compute ancestry proportions (physical or genetic) pA and pB given a list of 4 haplodicts. First 2 represent the parent corresponding to larger t , while the latter two correpsond to parent with smaller t. 
    '''
#{{{ 
    assert len(haplodicts) == 4

    pAs,pAs_1 = [],[]
    for hap_dicts in haplodicts:
        assert set([x['haplotypes'] for x in hap_dicts]) <= set(['0','1'])
        total = sum([x['length'] for x in hap_dicts])
        l0 = sum([x['length'] for x in hap_dicts if x['haplotypes'] == '0'])
        l1 = sum([x['length'] for x in hap_dicts if x['haplotypes'] == '1'])
        pAs.append(l0/total)
        pAs_1.append(l1/total)
 #}}}
    return(pAs,pAs_1)
        
#### OPTIONS #####
def collapse_blocks(haplodicts, alternate_key=None):
    '''
    Collapses contiguous block in haplo-parsed dictionary
    '''
#{{{

    if alternate_key:
        range_key=alternate_key
    else:
        range_key='startstop'

    #Loop over blocks in haplodicts
    contig_list = []
    extending = False
    for index,elem in enumerate(haplodicts):  
    
        #While extending, if next haplotype is not the same as current OR not on the same chromosome, end extended block at current, else go to next block
        if extending:

            #First we check if we've reached the end 
            if index == len(haplodicts) - 1: #If we've extended into last block, save end state and continue
                contig_block_end = elem[range_key][1]

                new_elem = {'haplotypes':contig_block_haplotype,range_key:[contig_block_start,contig_block_end], 'chromosome':contig_block_chromosome, 'length':contig_block_end - contig_block_start, 'end_status':elem['end_status']}

                contig_list.append(new_elem)

                extending = False
                continue

            if elem['haplotypes'] != haplodicts[index+1]['haplotypes'] or elem['chromosome'] != haplodicts[index+1]['chromosome'] : #If we've extended into different block (different chromosome or different haplotype ), save end state and continue 
                contig_block_end = elem[range_key][1]

                new_elem = {'haplotypes':contig_block_haplotype,range_key:[contig_block_start,contig_block_end],'chromosome':contig_block_chromosome,'length': contig_block_end - contig_block_start , 'end_status':elem['end_status']}

                contig_list.append(new_elem)

                extending = False
                continue

            else:
                continue

        #If we've reached the end, and are not extending, append current element to final list 
        if index == len(haplodicts) - 1:
            contig_list.append(elem)
            continue

        #If next block is same as current, AND on the same chromosome :  Save start, chromosome and haplotype information and start extending. Otherwise append current element to final list and continue  
        if elem['haplotypes'] == haplodicts[index+1]['haplotypes'] and elem['chromosome'] == haplodicts[index+1]['chromosome']:     
            contig_block_start = elem[range_key][0]
            contig_block_haplotype = elem['haplotypes'] 
            contig_block_chromosome = elem['chromosome'] 
            extending = True
        else:
            contig_list.append(elem)
            continue
            #}}}
    return(contig_list)

