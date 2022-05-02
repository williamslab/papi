import re

def splitbp_lines(bplines, outfile_prefix):
    '''
    Args:
    bpfile : Input admixsimu format bpfile  ( or gen file)

    Returns:
    NULL. Writes out admixSimu format bp/gen files ( one for each chromosome ) 
    '''

#{{{     
    #Loop over chromosomes
    for chrom in range(1,23):
        regex_chr = '\s'  + str(chrom) + '\|\s((\d+:\d+(\s|$))+)'
        writelines=[]
        for bpline in bplines:
            match = str(re.search(regex_chr,bpline).group(1))[:-1]
            writelines.append(match) #Get 1st capturing group, remove trailing whitespace 
        outfile = outfile_prefix[:-7] + '_chr{}.asu.bp'.format(chrom)
        with open(outfile,"w") as writefile:
            writefile.write("YRI CEU\n")
            writefile.write("\n".join(map(str,writelines)))
#}}}
    return(writelines)


if __name__ == '__main__':

    import argparse
    import re

    parser = argparse.ArgumentParser(description="takes in a asu.gen(bp) file and splits it into chromosomes, giving us the final admixsimu formatted files")

    # add arguments
    parser.add_argument("--inputfile", "-i", required=True, help="set input file (must be a bp file)")

    # read arguments from the command line
    args = parser.parse_args()
    # read static arguments 

    with open(args.inputfile,'r') as f:
        inlines = f.readlines()
        
    splitbp_lines(inlines, args.inputfile) 
