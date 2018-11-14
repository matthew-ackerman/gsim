POPULATION=4800				# The total population size being simulated
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
REF="reference.fa"  			# Name of the reference from which reads are simulated, this can be anything you want.
ASSEMBLY="assembly.fa"			# Name of the assembly to which reads are aligned. 
COV=15					# Average sequencing coverage (Poisson distributed).
SNPS=10000				# Number of loci to simulate in pedigree_sim. This doesn't have to be the same as the actual reference genome size. 
VAR=320000
REFTYPE='D'				# Reference types are 'Y'east, 'D'rosophila, and 'R'andom

LD_DIST=1000				#

STRUCTURE="g"				# for possible structures see pedigree_sim -h

F=1					# inbreeding value

EA=0					# Expectation of additive genotypic effects.
ED=0					# Expectation of dominance genotypic effects.
VA=1					# Variation of additive genotypic effects
VD=0					# Variance of dominance genotypic effects
VAD=0					# Covarinace between additive and dominance genotypic effects.
VE=0.5					# Variance of environmental effects.
