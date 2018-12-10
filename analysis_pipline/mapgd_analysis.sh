#!/bin/bash

cd ..
source settings.sh
cd /nfs/users/nfs_m/ma18/src/gsim/analysis_pipline

samtools view -H ../sequences/seq_000${paired}${filtered}${realign}${clipped}${bwasuffix} > ../sequences/temp-header.txt

name=$1
LD_DIST=$2

samtools mpileup -B ../sequences/*${paired}${filtered}${realign}${clipped}${bwasuffix} -f $name | gzip - > ../sequences/mpileup.txt.gz
mapgd proview -H ../sequences/temp-header.txt -n ../sequences/name-file.txt -s | gzip - > ../sequences/pro.txt.gz
echo "calling alleles"
mapgd allele -i ../sequences/pro.txt.gz -c 1 -g 2 -e 0.0001 -b | mapgd filter -q 0.001 -p 6 -g 2 -c 300 -C 3000 | gzip - > ../sequences/mapgd_calls.txt.gz
zcat ./mapgd_calls.txt.gz | mapgd filter -q 0.001 -p 6 -g 2 -X 0.05 | gzip - > ../sequences/mapgd_calls_p6_g2_X_0.05.txt.gz
exit


#mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 20 -e 0.0001 | mapgd filter -q 0.001 -p 1 -g 2 -N 1 | gzip - > ../analysis_files/mapgd_calls.txt.gz
zcat ../analysis_files/mapgd_calls.txt.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_calls-trim.csv
mapgd allele -i ../analysis_files/pro.txt.gz -c 1 -g 2 -e 0.0001 | mapgd filter -q 0.01 -p 10 -g 10 -N 1 | gzip - > ../analysis_files/mapgd_calls.txt.gz
echo "estimating genotypes"
mapgd genotype -p ../analysis_files/pro.txt.gz -m ../analysis_files/mapgd_calls.txt.gz | gzip - > ../analysis_files/genotype.gcf.gz
echo "estimating ld"
mapgd linkage -i ../analysis_files/genotype.gcf.gz -D $LD_DIST | gzip - > ../analysis_files/mapgd_linkage.out.gz
zcat ../analysis_files/mapgd_linkage.out.gz | tail -n +6 | sed '$d' >  ../analysis_files/mapgd_linkage-trim.csv
echo "estimating relatedness"
cat ../analysis_files/genotype.gcf.gz | gunzip - | mapgd relatedness > ../analysis_files/mapgd_relatedness.out
echo "estimating quantitive componenets"

#mapgd quant -r ../analysis_files/mapge_relatedness.out -p ../pedigree.txt 
