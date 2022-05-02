import pandas as pd 
import numpy as np 
import os 
import random

def distribute(oranges, plates):
    base, extra = divmod(oranges, plates)
    return [base + (i < extra) for i in range(plates)]


if __name__ == '__main__':

#{{{ 

    import argparse
    
    parser = argparse.ArgumentParser()
    
    # add arguments
    parser.add_argument("--inputdir", "-I", required=True, help="set input directory containing .asu.mrk.bp files")
    parser.add_argument("--t1", "-t1",type=int, required=True, help="t1")
    parser.add_argument("--t2", "-t2",type=int, required=True, help="t1")
    parser.add_argument("--splits", "-s", type=int, required=True,help="specifies no of splits as well as split ratio (e.g splits=4 will give 4 sets of splits with 1:3 simulation set to hold set as sim.txt and hold.txt files )")
    parser.add_argument("--ninds", "-n", type=int,required=True,help="number of individuals in each split")
    parser.add_argument("--indsfile","-i", required=True, help="inds file output from generate_admixsimuBpfiles.sh")
    parser.add_argument("--outdir","-o", required=True, help="Outfile prefix")

    #read arguments from the command line
    args = parser.parse_args()


    indsDF_full = pd.read_table(args.indsfile,header = None)  
    #Read focal entries from indsDF_full
    indsDF_foc = indsDF_full

    total_inds = args.ninds*args.splits
    assert len(indsDF_foc) == total_inds*2

    #Split indsDF_foc into args.n splits
    indsDF_final_splits = np.array_split(indsDF_foc,args.splits)     
    
    #Write out the n splits into .txt files
    for i,df in enumerate(indsDF_final_splits):
        df2 = df[df.index % 2 == 0] #We only want every other line for writing
        indfilepath,indfilename = os.path.split(args.indsfile)

        writefilename = args.outdir + '/' + indfilename[:-5] + '_split{}'.format(i+1)  + '.inds'
        df2.to_csv(writefilename,sep="\t", header=False,index=False)

    #Subset each chromosomal .asu.mrk.bp file into the splits using splits_indices 
    for chrom in range(1,23):

        infile = args.inputdir + '/{}/'.format(chrom) + str(args.t1) + '_' + str(args.t2) +'gens_chr{}'.format(chrom)+ '.asu.mrk.bp'
        with open(infile,'r') as f:
            inlines = f.readlines()
            
        inlines_wo_header = inlines[1:]
        header = inlines[0]

        for split in [1,2,3,4]:
            outlines = [header] + inlines_wo_header[(split-1)*44:split*44] 
            path, infilename = os.path.split(infile) 
            outfile = args.outdir + '/{}/'.format(chrom) + infilename[:-11] + '_splits{}'.format(split) + '.asu.mrk.bp' 
            with open(outfile,'w+') as f:
                f.writelines(outlines)
#}}}

