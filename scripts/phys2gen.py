import re   
import pandas as pd 


def phys2gen_file(bpfile,mapfile):

    '''
    Args:
    bpfile : Input bpfile in standard format with physical positions
    mapfile : mapfile with mapping from physical position to genetic position

    Returns:
    NULL. Writes a new bpfile (bpfile.gen) with genetic instead of physical positions
    '''

#{{{1

    #Read in mapfile as DataFrame 
    mapDF = pd.read_table(mapfile) 
    outlist = []

    #Loop over bp file line by line 
    with open(bpfile) as f:
#{{{2 
        for count, bpline in enumerate(f):
            bpline = bpline.replace('\n' , ' ') #Replace trailing newline with space for easire regex processing

            parsed_list = re.findall(r'(\d{1,2}\|\d+\s(?:[A-Za-z0-9]+:\d+\s)+)', bpline, re.MULTILINE) #break up string into chromosomes
            header = re.findall(r'(.*\s)1\|', bpline, re.MULTILINE) #Store header to write to output later

            #Loop over chromosomes, get physical positions, interpolate and replace 
            write_list = [None]*len(parsed_list)
            for i in range(len(parsed_list)):
#{{{ 
                chrom = re.search(r'^[^|]*',parsed_list[i]).group(0)
                chrom_start = int(re.search(r'\|(.*?)\s', parsed_list[i]).group(1))
                chrom_end_list = re.findall(r':(\d+)', parsed_list[i])

                #make chrom_start and chrom_end_list into one list and convert to ints 
                physical_positions = [chrom_start] + [int(x) for x in chrom_end_list]

                #Subset dataframe ( current chromosome ) 
                mapDF_chrom = mapDF[mapDF['#chr'] == int(chrom)].reset_index()


                #loop over physical positions 
                for p in physical_positions:

                    #Binary search DataFrame and find positions to interpolate between 
                    indx = mapDF_chrom.pos.searchsorted(p)
                   
                    #Do the interpolation
                    if not p == int(mapDF_chrom.pos[indx]): #if position is not found exactly 
                        num = float(mapDF_chrom.sexavg[indx]) - float(mapDF_chrom.sexavg[indx-1])
                        denom = int(mapDF_chrom.pos[indx]) - int(mapDF_chrom.pos[indx-1])
                        slope = num/denom
                        geneticpos_interp = round(float(mapDF_chrom.sexavg[indx-1]) + slope*(p - int(mapDF_chrom.pos[indx-1])),3)
                     
                    else:
                        geneticpos_interp = round(float(mapDF_chrom.sexavg[indx]),3)  
                    #Replace current physical position with genetic position in parsed_list[i] 
                     
                    parsed_list[i] = parsed_list[i].replace('|' + str(p), '|' + str(geneticpos_interp))              
                    parsed_list[i] = parsed_list[i].replace(':' + str(p), ':' + str(geneticpos_interp))              

                hap_strs =  re.findall(r'\s\d+:', parsed_list[i]) 
                genpos_strs =  re.findall(r':\d+\.\d+\s', parsed_list[i]) 

                #Check that output makes sense
                assert all(['.' not in x[2:-1] for x in hap_strs])
                assert all(['.' in x[2:-1] for x in genpos_strs])
#}}}

            #Smoosh the lists of strings together with the header and add newline
            joinlist = header + parsed_list 
            outlist.append(''.join(joinlist)) 
#}}}2


    #Write out modified bpfile 
    outfile = bpfile[:-3] + '.gen'
    with open(outfile, 'w') as f:
        f.write("\n".join(map(str,outlist)))

#}}}1

    return()               


if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser()

    # add long and short argument
    parser.add_argument("--inputfile", "-i",required=True, help="set input file (must be a bp file)")
    parser.add_argument("--mapfile", "-m",required=True, help="set mapfile path")
    # read arguments from the command line
    args = parser.parse_args()

    res = phys2gen_file(args.inputfile, args.mapfile)

