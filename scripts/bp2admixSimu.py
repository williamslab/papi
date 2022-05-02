if __name__ == '__main__':

    import argparse
    import pandas as pd 
    import numpy as np 
    import admixString_generators as strGen
    import parsers as prs
    import re
    import random

    parser = argparse.ArgumentParser(description="takes in a bp file and creates admixsimu format bp files with ancestry states for all possible ancestry strings and t1,t2values consistent with the input bpgen file. Can be run with or without subadmixture.Also creates a .txt table file with bp-headers, admixture strings, p1, p2, t1,t2 columns.")
    # add arguments
    parser.add_argument("--inputfile", "-i", required=True, help="set input file - must be a bp file")
    parser.add_argument("--t1", "-t1", required=True, type=int, help="number of observed meioses on the larger t branch")
    parser.add_argument("--t2", "-t2", required=True, type=int,help="number of observed meioses on the smaller t branch")
    parser.add_argument('-o','--options', nargs='+', help='options for subadmixture - e.g -o n s - for normal and subadmixed on the respective branches', required=True)

    # read arguments from the command line
    args = parser.parse_args()
    # read static arguments 
    
#{{{ 

    #Generate 4 sets of admixstrings ( one for each split : 50/50, 40/60, 25/50,  25/75 )

    admixstr_list = []
    for ptup in [(0.2,0.2),(0.2,0.4),(0.2,0.6),(0.2,0.8)]:
        admixstr_list.append(strGen.generateAdmixStrings(args.t1, args.t2,ptup[0],ptup[1], options=[args.options[0],args.options[1]]))
    
    print('Generated strings!')
    with open(args.inputfile, 'r') as f:
        lines_all = f.readlines()
        
    nfams= (len(lines_all)//6) #Batch size of individuals for each admixtype. 2 haplotypes for the individual, and 4 for his/her parents

    tablelines = [['header\tadmixstring\tpA\tpB\tt1\tt2']]
    writelines = []

    for fam in range(nfams):
        
        admixstrs = admixstr_list[fam//22]
        admixstr = admixstrs[fam%len(admixstrs)] 
        start = fam*6
        end = start + 6
        famlines = lines_all[start:end]

        famlines_anc = [] 
        for j in range(len(famlines)):

            if j%2 != 0: 
                continue
            famlines_anc.append(prs.ancestryParse([famlines[j], famlines[j+1]],list(admixstr)))

        famlines_anc_flat = [x for y in famlines_anc for x in y]
        famlines_haplodicts = [prs.haplo_parse(x) for x in famlines_anc_flat]

        print('famlines haplo parsed')

        ps_0,ps_1 = prs.computeAncprops(famlines_haplodicts[:4])
        assert 0.99 <= ps_0[0] + ps_1[0] <= 1.01, "t1 = {} , t2 = {}, admixstr = {}, header={}".format(str(args.t1),str(args.t2), admixstr, famlines_anc_flat[0].split(' ')[0])

        pA = sum(ps_0[:2])/2
        pB = sum(ps_0[2:])/2

        tablelines.append([x.split(' ')[0] +'\t'+ admixstr + '\t' +str(pA) + '\t' + str(pB) + '\t' + str(args.t1) + '\t' + str(args.t2)  for x in famlines_anc_flat])
        writelines.append(prs.admixSimu_parse(famlines_anc_flat[4:]))


    tablefile = str(args.t1) + '_' + str(args.t2) + '_' + 'pois_sxavg.inds'

    with open(tablefile, 'w+') as f:
        f.write("\n".join([x for y in tablelines for x in y]))

    writelines_flat = [x for y in writelines for x in y]
    outfile = str(args.t1) + '_' + str(args.t2) + '_' + 'pois_sxavg.asu.bp'
    with open(outfile, 'w+') as f:
        f.write("".join(writelines_flat))

        #Write file for current typ
#        infile = re.search('\/([^\/]*$)', args.inputfile).group(1)
#}}}
