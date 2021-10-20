import ast 
import sys 

def read_data(filename, col_num):
    '''
    Read in the chosen column number of input file and literal_eval it. Outputs a list of literal_eval'ed entries corresponding to chosen column  
    '''
    datalist = [] 

    with open(filename) as f:
        for line in f:
            entry_string = line.split('\t')[col_num - 1]
            if entry_string[0] != 't':
                entry = ast.literal_eval(entry_string)
            else:
                continue
            datalist.append(entry)

    return(datalist)

def read_data_raw(filename, cols):
    '''
    Read in the chosen column numbers of input file. Outputs a list strings ( with tab separators) corresponding to chosen columns  
    '''
    datalist = [] 

    with open(filename) as f:
        for line in f:

            line_list = line.split('\t')
            entry_list = [line_list[x-1] for x in cols]

            if entry_list[0] != 'header':
                entry_string = '\t'.join(entry_list)
            else:
                continue

            datalist.append(entry_string)

    return(datalist)

def parseUnordered_dataDicts(dataDicts):
    '''
    Parse list of lists of dicts containing ordered tract information in chromosomes. Output is dict of form {'11':304.5, '10':103.4, '00':212.2 } containing tract type and total morgan length of genome in tract of that type.
    '''
    result = {'11':0, '10':0, '00':0}

    for elem in dataDicts:

        if elem['type']=='11':
            result['11'] += elem['length']
        elif (elem['type']=='10') | (elem['type']=='01'):
            result['10'] += elem['length']
        elif elem['type']=='00':
            result['00'] += elem['length']
        else:
            sys.exit('Couldnt find appropriate tract type')
     
    
    #Check that all keys are present 
    assert sorted(list(result.keys())) == sorted(['00','11','10'])


    return(result)

def parseClean_dataDicts(dataDicts):
    """
    Outputs more sensible list of dicts representation of tract information. Used for evaluating pooled markov process likelihood. Input is list of lists of dicts representing tract information arranged in chromosomes
    """
    dataDicts_clean = [{'10':dicT['01']} if list(dicT.keys())[0] == '01' else dicT for dicT in dataDicts]#Remove redundant keys


    return([{'type':list(x.keys())[0], 'length':list(x.values())[0]} for x in dataDicts_clean])

