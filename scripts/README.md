# Prerequisites

You must have [ped-sim](https://github.com/williamslab/ped-sim) and [admix-simu](https://github.com/williamslab/admix-simu) installed. Additionally a set of phased YRI and CEU genotypes (at least 88 of each) in eigenstrat format to use both as admix-simu panels and for local ancestry inference with HapMix.

The [convertf](https://reich.hms.harvard.edu/software/InputFileFormats) utility can be used to convert plink format data to eigenstrat format.


# Running the scripts 

## 1. generateAdmixsimuFiles.sh

```generateAdmixsimuFiles.sh``` takes in 4 user supplied arguments which are t1,t2, a genetic mapfile, a panel-directory containing eignestrat format phased panel files(see above), and an output directory for results.  Note that the panel-directory must contain a set of per-chromosome .snp files(eigenstrat format) must be supplied by the user and named YRIpanelfile_${chr}.snp. Genetic distance information is not strictly required in these snp files, so the script should run even in its absence.

Here's an example usage string  -
```
sh generateAdmixsimuFiles.sh 9 5 hapmap_mapfile.txt panelDir results
``` 
Which will generate output for t1=9, t2=5 using the ```hapmap_mapfile.txt``` mapfile and will store the output in ```results```.

## 2. admixsimu_simulate.sh

```admixsimu_simulate.sh``` takes in 4 user supplied arguments which are t1,t2, an output directory, and a directory that contains eigenstrat format phased panel files.

This creates the appropriate eigenstrat files - .geno for unphased genetic information and a .inds and .rates file required by hapmix (The snp2rates.py script requires a snp file with genetic information - so ensure that this column isn't empty ).



