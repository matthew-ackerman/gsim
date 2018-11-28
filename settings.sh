###### Simulation Settings


POPULATION=1200				# The total population size being simulated
POP2=$((POPULATION/2))			# Population size being simulated after a reduction in population size
# 
# Note: whether values like POP2 have any impact on the simulation will be determined by the kind of simulation that is being run, which in tern depends on parameters passed to pedigree_sim. 
# The most important value passed to pedigree_sim is the argument passed to pedgigree_sim --type.
#
SAMPLE=75				# Number of individuals sampled from the population for sequencing.
SAMP=`seq 0 1 $((POPULATION-1)) | shuf | head -$SAMPLE` # Simulated sample names. This is just a list from 0 to population size, and you probably shouldn't change that 

STATE_FILE="`readlink -f ./sequences`/state.txt"

for s in $SAMP
do
SAMPLE_CHRM+=","
SAMPLE_CHRM+=$((s*2+2))
SAMPLE_CHRM+=","
SAMPLE_CHRM+=$((s*2+3))
SAMPLE_NAME+=","
SAMPLE_NAME+=$((s+2))
done

TIME=$((5*POPULATION))			# Number of generation simulated.
TIMEX=$((9*PULATION+POPULATION/2))	# Time of the change between POPULATION and POP2
REF=`readlink -f ./sequences`"/reference.fa"  			# Name of the reference from which reads are simulated, this can be anything you want.
ASSEMBLY=`readlink -f ./sequences`"/assembly.fa"			# Name of the assembly to which reads are aligned. 
ASSEMBLYID="001"
COV=8					# Average sequencing coverage (Poisson distributed).
SNPS=20400				# Number of loci to simulate in pedigree_sim. This doesn't have to be the same as the actual reference genome size. 
VAR=300
REFTYPE='D'				# Reference types are 'Y'east, 'D'rosophila, and 'R'andom

LD_DIST=1000				#

STRUCTURE="g"				# for possible structures see pedigree_sim -h

F=32					# inbreeding value

EA=0					# Expectation of additive genotypic effects.
ED=0					# Expectation of dominance genotypic effects.
VA=1					# Variation of additive genotypic effects
VD=0					# Variance of dominance genotypic effects
VAD=0					# Covarinace between additive and dominance genotypic effects.
VE=0.5					# Variance of environmental effects.

###### Variant Calling Settings. These need to be shared across all variant calling pipelines.

Q=7					# Read based filtering flags to be passed to samtools, if applicable.
q=7					# (Mapping quality)

WD=`readlink -f ./sequences/`

###### Programs 

picard="java -jar "`readlink -f ../picard/picard.jar`
Trimmomatic="java -jar "`readlink -f ../Trimmomatic-0.38/trimmomatic-0.38.jar`
GATK="java -jar "`readlink -f ../gatk/GenomeAnalysisTK.jar`
bamUtil=`readlink -f ../bamUtil/bin/bam`
bwamem="bwa mem"


##### 

assemblyDir=$WD
assemblyName=$ASSEMBLY

fastqBase=$WD
SampleDir="" #$WD
CloneID=KAP-00074

#this should be different from ? if you want to simulate removing the wrong adaptors.
adapterFile=`readlink -f ./sequences/adapters.fa`
assemblyID=PA42_4_0

GENUS="Cacaradon"

suffix=".fq.gz"
forward=".1"
reverse=".2"
paired=".paired"
unpaired=".unpaired"
sorted=".sorted"
filtered=".filtered"
realign=".realigned"
clipped=".clipped"
bwasuffix=".bam"

sequenceDict=${assemblyName%.fa}.dict
