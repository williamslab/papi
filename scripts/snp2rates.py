import pandas as pd 
import numpy as np

def snp2rates(snpfile):

    '''
    converts snp file to rates file

    Usage: python ~/siddharth/Tools/pythonpackages_sid/snp2rates.py {vcffile}
    '''
#{{{ 
    with open(snpfile) as f:
        snpfile_lines = f.readlines()

    physical_positions = [row.split()[3] for row in snpfile_lines ]
    genetic_positions = [row.split()[2] for row in snpfile_lines ]
    nsites = len(physical_positions)

    ratesfile = snpfile[:-3] + 'rates'

    with open(ratesfile,'w') as f:
         f.write(':sites:' + str(nsites) + '\n' + ' '.join(physical_positions) +'\n' + ' '.join(genetic_positions))
#}}}

    return()

if __name__ == "__main__":

    import sys

    snp2rates(sys.argv[1])


