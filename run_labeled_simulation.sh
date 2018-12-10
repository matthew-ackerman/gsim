#!/bin/bash

source settings.sh

PED_SIM=`readlink -f ./pedigree_simulation/pedigree_sim`

:'
echo "simulating population"

cd sequences

echo $SAMPLE
echo $POPULATION

$PED_SIM -y $STRUCTURE -f $F -N $POPULATION -1 $SAMPLE -s $SNPS -g $TIME -plnGt -k $TIME -v $VAR -e 0.005 -t > states.txt 

#head -2 name-file.txt > temp-file.txt
#head -3 name-file.txt | tail -n 1 | cut -d '	' -f 1,$SAMPLE_NAME >> temp-file.txt
#tail -n 1 name-file.txt >> temp-file.txt
#mv temp-file.txt name-file.txt

rm -rf pedigree.txt.gz
	gzip pedigree.txt
cd ..

case $REFTYPE in
  [S]   )
		echo "simulating reference"
		python reference_simulation/mutation_simulation2.py -l 24 -S $SIZE > ./sequences/reference_mutations.txt
		python reference_simulation/mutation_simulation.py -m ./sequences/reference_mutations.txt -s ./reference_simulation/seed.fa > $REF
		;;
  [Y]   )
		echo "using yeast chromosome IV"
		zcat ./real_genomes/S288C_Chromosome\ IV.fsa.gz > $REF
		zcat ./real_genomes/S288C_Chromosome\ IV.fsa.gz > $ASSEMBLY
                zcat ./real_genomes/TruSeq.fa.gz > $adapterFile
		;;
  [D]   )
                echo "using dmel3R "
                zcat ./real_genomes/dmel-3R.fa.gz > $REF
                zcat ./real_genomes/dmel-3R_2003.fa.gz > $ASSEMBLY
#                zcat ./real_genomes/dmel-3R.fa.gz > $ASSEMBLY
                zcat ./real_genomes/TruSeq.fa.gz > $adapterFile
                ;;
esac


echo "making individual genomes"
cd variant_simulation
bash state_to_fasta.sh $REF ../sequences/states.txt
cd ..
echo "simulating sequencing"
:'
#rm -rf ./sequences/states.txt.gz
#gzip ./sequences/states.txt

for x in $(seq -f "%03g" 0 1 $((SAMPLE-1)) )
do 
	NAME=seq_$x

	echo $NAME

	gunzip ./sequences/$NAME.0.fa.gz
	gunzip ./sequences/$NAME.1.fa.gz

	./sequencing_simulation/art_illumina --id "0.fa_" -qs -15 -qs2 -15 -ss HS25 -sam -i ./sequences/$NAME.0.fa -p -l 150 -f $COV -m 200 -s 10 -o ./sequences/temp.0 > /dev/null
	./sequencing_simulation/art_illumina --id "1.fa_" -qs -15 -qs2 -15 -ss HS25 -sam -i ./sequences/$NAME.1.fa -p -l 150 -f $COV -m 200 -s 10 -o ./sequences/temp.1 > /dev/null

	gzip ./sequences/$NAME.0.fa
	gzip ./sequences/$NAME.1.fa

	cat ./sequences/temp.01.fq | gzip - > ./sequences/$NAME${forward}.fq.gz
	cat ./sequences/temp.11.fq | gzip - >> ./sequences/$NAME${forward}.fq.gz
	cat ./sequences/temp.02.fq | gzip - > ./sequences/$NAME${reverse}.fq.gz
	cat ./sequences/temp.12.fq | gzip - >> ./sequences/$NAME${reverse}.fq.gz

	rm ./sequences/temp.01.fq
	rm ./sequences/temp.11.fq
	rm ./sequences/temp.02.fq
	rm ./sequences/temp.12.fq


	rm ./sequences/temp.01.aln
	rm ./sequences/temp.11.aln
	rm ./sequences/temp.02.aln
	rm ./sequences/temp.12.aln

	samtools sort ./sequences/temp.0.sam -b  > ./sequences/$NAME${forward}.true.bam
	samtools sort ./sequences/temp.1.sam -b  > ./sequences/$NAME${reverse}.true.bam

	rm ./sequences/temp.0.sam
	rm ./sequences/temp.1.sam

done 

cd variant_calling_pipeline
bash bwa_pipeline.sh 
cd .. 

ls sequences/ | grep  ${paired}${filtered}${realign}${clipped}${bwasuffix} > sequences/bam_list.txt

#mv ./sequenes/polymorphism.map
#python injection.py -a ./sequences/merged_seq.bam -r ./sequences/merged.sort.bam -s ./sequences/polymorphisms_ref.map > ./sequences/polymorphisms.map


cd analysis_pipeline

#./mapgd_analysis_newton.sh $REF
./mapgd_analysis.sh $ASSEMBLY $LD_DIST
exit
./bcftools_analysis.sh $ASSEMBLY
./angsd_analysis.sh $ASSEMBLY
./gatk_analysis.sh $ASSEMBLY
./plink_analysis.sh $ASSEMBLY
./gcta_analysis.sh 

gunzip ../sequences/states.txt
python get_frequencies.py ../sequences/states.txt ../sequences/polymorphisms.map > ../analysis_files/true_frequencies.csv
#python get_ld.py ../sequences/states.txt ../sequences/polymorphisms.map $LD_DIST > ../analysis_files/true_ld.csv
gzip ../sequences/states.txt

Rscript Ackerman2017/make_figure_1a.rscript	#Bias RMSE of allele frequencies
Rscript Ackerman2017/make_figure_1c.rscript	#Bias RMSE of inbreeding
Rscript Ackerman2017/make_figure_1d.rscript	#Bias RMSE of LD.
Rscript Ackerman2017/make_figure_1e.rscript	#ROC.

#./mapgd_benchmark.sh $ASSEMBLY > ../analysis_files/benchmark.csv
#./bcftools_benchmark.sh $ASSEMBLY >> ../analysis_files/benchmark.csv
#./angsd_benchmark.sh $ASSEMBLY >> ../analysis_files/benchmark.csv
#./gatk_benchmark.sh $ASSEMBLY >> ../analysis_files/benchmark.csv
