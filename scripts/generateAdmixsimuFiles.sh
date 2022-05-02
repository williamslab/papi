#!/bin/bash

t1=$1
t2=$2
MAP_FILE=$3 #mapfile for ped-sim run 
OUT_DIR=$4 #specify output directory for all files

##### ped-sim run #######
DEF_FILE=${t1}gens.def
t=$((t1+1)) #Definition of t in pedsim is different to our definition of number of observed crossovers
python ./deffileGenerator.py ${t} > ${DEF_FILE} #generates ped-sim format definition file

echo "${t1} gens is being simulated with def file:${DEF_FILE}"
ped-sim --bp -d ${DEF_FILE} -m ${MAP_FILE} -i /dev/null --pois -o ${OUT_DIR}/${t1}gens.bp --founder_ids 


#python src/d00_data/phys2gen.py -i ${t1}pois_sxavg.bp -m ${MAP_FILE} #Run this to obtain bp files with genetic instead of physical distance


######## Creating admixsimu format files #########
echo "Runninng bp2admixSimu for t1 = ${t1} and t2 = ${t2}"
python src/d00_data/bp2admixSimu.py  -i ${OUT_DIR}/${t1}gens.bp -t1 ${t1} -t2 ${t2} -o n n #creates files with extension .asu.bp


####### Split files into individual chromsome files #######
for file in ${OUT_DIR}/*.asu.bp
do
    echo ${file}
    python ./splitChrom_AdmixSimu.py -i ${file}
done

######## Reformat .bp files with marker position instead of physical postiion
#WARNING - This step assumes that the user has eigenstrat format .snp files derived from the panel data (see README)
for chr in {1..22} 
do
    python ./phys2mark_AdmixSimu.py -i ${OUT_DIR}/${t1}_${t2}gens_chr${chr}.asu.bp -s snpfile${chr}.txt
done

####### Create splits ###############
gen=$(($t1+2))
grep "_g$gen-" ${OUT_DIR}/${t1}_${t2}gens.inds > ${OUT_DIR}/${t1}_${t2}gens_foc.inds
python ./splits_AdmixSimu.py -I ${OUT_DIR} -t1 "$t1" -t2 "$t2" -s 4 -n 22 -i ${OUT_DIR}/${t1}_${t2}gens_foc.inds -o ${OUT_DIR}


