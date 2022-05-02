import re
import numpy as np
import pandas as pd 
def phys2mark_lines(inlines,snpDF):

    snpDF.columns = ['rsid','chr','sexavg','pos','ref','alt']
    snpDF_pos = np.array(snpDF["pos"])
    outlines=[]
    for inline in inlines:
        haplist = re.findall(r'(\d):', inline)
        physlist = [int(x) for x in re.findall(r':(\d+)', inline)]
        markers = [x + 1 for x in np.searchsorted(snpDF_pos, physlist)]
        outline = " ".join(["{}:{}".format(x[0],x[1]) for x in list(zip(haplist,markers))]) 
        outlines.append(outline)

    return(outlines)

if __name__ == '__main__':

    import argparse
    import re

    parser = argparse.ArgumentParser(description="takes in a single chromosome asu.bp file and replaces physical position with marker number everywhere")

    # add arguments
    parser.add_argument("--inputfile", "-i", required=True, help="set input file (must be a bp file)")
    parser.add_argument("--snpfile", "-s", required=True, help="set snpfile")
    # read arguments from the command line
    args = parser.parse_args()
    # read static arguments 

    with open(args.inputfile,'r') as f:
        inlines = f.readlines()
        
    if inlines[0] == "YRI CEU\n":
        header = inlines[0]
        inlines=inlines[1:]

    snpDF = pd.read_table(args.snpfile)

    outlines = phys2mark_lines(inlines, snpDF)
    writefile = args.inputfile[:-2] + 'mrk.bp'
    with open(writefile,'w+') as f:
        f.write(header)
        f.writelines("\n".join(outlines))


