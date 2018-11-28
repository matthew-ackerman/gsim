######### Loading carbonate-specific resources (change if not using this cluster) #######

cd ..
source settings.sh
cd /nfs/users/nfs_m/ma18/src/gsim/variant_calling_pipeline/

# please install ngsutils (see INSTALL.txt for instructions)
# then, provide the path to the fastqutils binary 
# please note that fastqutils is only required for the original pipeline (i.e. the one using novoalign)
# fastqutils=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall/software/ngsutils/bin/fastqutils

####### Path to adapters #######

###### Misc names ######


##### Number of threads #####
nThreads=8

#a="

cd $WD


echo "Trimming adapter sequences from sequence reads."

#This command needs to make an array of all the root names for reads.
NAMES=(`ls | grep -o 'seq_[0-9]\{3\}.1' | grep -o 'seq_[0-9]\{3\}'`)

echo "Making a dictionary file of a reference."
$picard CreateSequenceDictionary R=$ASSEMBLY O=$sequenceDict

for i in "${!NAMES[@]}" 
do
seq_f=${NAMES[${i}]}${forward}
seq_r=${NAMES[${i}]}${reverse}
${Trimmomatic} PE ${seq_f}${suffix} ${seq_r}${suffix} ${seq_f}${paired}${suffix} ${seq_f}${unpaired}${suffix} ${seq_r}${paired}${suffix} ${seq_r}${unpaired}${suffix} HEADCROP:3 ILLUMINACLIP:$adapterFile:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30
done


samtools faidx $ASSEMBLY
bwa index $ASSEMBLY

for seq in ${NAMES[@]}
do
echo "$bwamem -t 16 -M $ASSEMBLY ${seq}${forward}${paired}${suffix} ${seq}${reverse}${paired}${suffix} | samtools view -b -q 20 -f 3 -F 3844 - " #> ${seq}${paired}${bwasuffix}
$bwamem -t 16 -M $ASSEMBLY ${seq}${forward}${paired}${suffix} ${seq}${reverse}${paired}${suffix} | samtools view -b -q 20 -f 3 -F 3844 - > ${seq}${paired}${bwasuffix}
#$bwamem -t 16 -M $ALIGNMENT ${seq}${forward}${unpaired}${suffix} | gzip - > ${seq}${foward}${unpaired}${suffix} 
#$bwamem -t 16 -M $ALIGNMENT ${seq}${reverse}${unpaired}${suffix} | gzip - > ${seq}${reverse}${unpaired}${suffix} 
echo "${seq}"
done

# 5. Sort the BAM file using Picard.
echo "Sorting the bam file using Picard."

for seq in ${NAMES[@]}
do
echo "$picard SortSam INPUT=${seq}${paired}${bwasuffix} OUTPUT=/dev/stdout SORT_ORDER=coordinate | $picard AddOrReplaceReadGroups INPUT=/dev/stdin OUTPUT=${seq}.temp${bwasuffix} RGID=${GENUS} RGLB=bar RGPL=illumina RGSM=${seq} RGPU=6"
$picard SortSam INPUT=${seq}${paired}${bwasuffix} OUTPUT=/dev/stdout SORT_ORDER=coordinate | $picard AddOrReplaceReadGroups INPUT=/dev/stdin OUTPUT=${seq}.temp${bwasuffix} RGID=${GENUS} RGLB=bar RGPL=illumina RGSM=${seq} RGPU=6
$picard MarkDuplicates READ_NAME_REGEX=null INPUT=${seq}.temp${bwasuffix} OUTPUT=${seq}${paired}${filtered}${bwasuffix} METRICS_FILE=${seq}.${ASSEMBLYID}_metrics.tx
$picard BuildBamIndex INPUT=${seq}${paired}${filtered}${bwasuffix}
done

# 9. Define intervals to target for the local realignment.

for seq in ${NAMES[@]}
do
$GATK -T RealignerTargetCreator -R ${ASSEMBLY} -I ${seq}${paired}${filtered}${bwasuffix} -o ${seq}.${assemblyID}.intervals
done

# 10. Locally realign reads around indels.


for seq in ${NAMES[@]}
do
$GATK -T IndelRealigner -R ${ASSEMBLY} -I ${seq}${paired}${filtered}${bwasuffix} -targetIntervals  ${seq}.${assemblyID}.intervals -o ${seq}${paired}${filtered}${realign}${bwasuffix}
done

# 11. Clip overlapping read pairs.
echo "Clipping the overlapping read pairs using bamUtil."

for seq in ${NAMES[@]}
do
$bamUtil clipOverlap --in ${seq}${paired}${filtered}${realign}${bwasuffix} --out ${seq}${paired}${filtered}${realign}${clipped}${bwasuffix} 
done

# 12. Index the clipped BAM file using Samtools
echo "Indexing the clipped BAM file using Samtools."
for seq in ${NAMES[@]}
do
samtools index ${seq}${paired}${filtered}${realign}${clipped}${bwasuffix} 
done
