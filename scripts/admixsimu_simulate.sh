#!/bin/bash

t1=$1
t2=$2
OUT_DIR=$3
paneldir=$4

for chr in {1..22}
do
    for split in {1,2,3,4}
    do
        echo "running admixsimu"
        echo ${t1},${t2},${split},${chr}
        
        #Run admixSimu
        ~/siddharth/programs/tools/admix-simu/simu-mix.pl -bp ${OUT_DIR}/${t1}_${t2}_gens_chr${chr}_splits${split}.asu.mrk.bp ${OUT_DIR}/${t1}_${t2}_pois_sxavg_chr${chr}_splits${split} -CEU ${paneldir}/CEUpanelfile_${chr}.phgeno -YRI ${paneldir}/YRIpanelfile_${chr}.phgeno  
       
        #Post processing admixsimu results to create HapMix format input files 
        echo "post-processing admixsimu results"
        python ./eig2eig.py ${t1}_${t2}_gens_chr${chr}_splits${split}.phgeno #create .geno file
        python ./snp2rates.py ${paneldir}/YRIpanelfile_${chr}.snp #create .rates file
        python ./geno2ind.py ${OUT_DIR}/${t1}_${t2}_pois_sxavg_chr${chr}_splits${split}.geno #Create .ind file
    done
done


