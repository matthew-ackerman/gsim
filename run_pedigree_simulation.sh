#!/bin/bash

source settings.sh

PED_SIM=`readlink -f ./pedigree_simulation/pedigree_sim`

echo "simulating population"

#64 is about right.

for F in 64 #8 16 32 64 128 256 512
do
cd sequences

$PED_SIM -y $STRUCTURE -f $F -N $POPULATION -1 $SAMPLE -s $SNPS -g $TIME -plnGtF -k $TIME -v $VAR -e 0.005 -t > $STATE_FILE

cd ..

mapgd writevcf2 -s $STATE_FILE -n ./sequences/name-file.txt | bcftools view -q 0.01 | gzip - > ./analysis_files/states_${F}.vcf.gz
mv ./sequences/cord-file.txt ./analysis_files/cord-file_${F}.txt
cat ./sequences/fstat.txt | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' > ./analysis_files/fstat_${F}.txt 
done
