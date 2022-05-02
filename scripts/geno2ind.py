def geno2ind(inputfile):

    '''
    Takes unphased genotype and creates a random ind file
    Write to stdout
    '''

#{{{ 
    with open(inputfile,'r') as f:
        inlines_one = f.readline()

    n_ind = len(inlines_one)

    writefile = inputfile[:-4] + 'ind' 
    with open(writefile,'w') as f:
        for i in range(1,n_ind):
            f.write('AA_{}\tM\tYC\n'.format(i))
#}}}



if __name__ == "__main__":

    import sys

    geno2ind(sys.argv[1])




