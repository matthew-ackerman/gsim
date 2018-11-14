#### Benchmark scripts.

angsd\_benchmark.sh
bcftools\_benchmark.sh
gatk\_benchmark.sh
mapgd\_benchmark.sh
	
These scripts should record the amount of real and system time a program spends analyzing simulated data with 1, 2, and 4 threads.

#### Analysis scripts.

angsd\_analysis.sh
bcftools\_analysis.sh
gatk\_analysis.sh
gcta\_analysis.sh
king\_analysis.sh
mapgd\_analysis\_newton.sh
mapgd\_analysis.sh
ldak\_analysis.sh
plink\_analysis.sh

breseq\_analysis.sh
mapgd\_pooled\_analysis.sh

ngsrelated.sh (Working ??)

mapgd\_analysis-vg.sh (Not working)
vg\_analysis.sh (Not working)

These scripts are called to estimate allele frequencies, f statistics, measures of pairwise relatedness, etc. from the simulated data.

#### Scripts that produce dot files.

draw\_pedigree\_focal\_indv.py
draw\_pedigree\_focal\_pair.py
draw\_pedigree.py

#### Heritability

gvcblup\_gwas.sh
gcta\_gwas\_setup.sh
gcta\_gwas.sh
plink\_gwas.sh

These scripts estimate the heritability of quantitative traits in the simulated data.

#### Subdirectories
mapgdr
pedigree
make\_figures
python\_utilities


#### Unclassified

closed\_form.py
count\_sub\_pedigree.py
draw\_benchmarks.py
estimate\_variance.rscript
factor\_2trait.py
factor5.py
factor.py
fast2.rscript
fast.rscript
figures.conf
get\_a\_d.py
get\_ad.py
get\_constants\_2.py
get\_constants.py
get\_frequencies\_from\_pol.py
get\_frequencies\_from\_vcf.py
get\_frequencies.py
get\_ld\_moments.py
get\_ld.py
get\_mu.py
get\_var.rscript
homozygosity\_track\_length.py
interpolate.py
lagrange\_factor2.py
lagrange\_factor.py
make\_mu.rscript
make\_qq\_from\_vcf.py
make\_qq\_plots.rscript
make\_r\_rel-fast-exp.rscript
make\_r\_rel-fast.rscript
make\_r\_rel.rscript
MKLMM.sh
mv\_sim.rscript
print\_last.py
py\_opt\_breeding.py
randfloat.py
rename\_markers.py
run\_labeled\_simulation.sh
states\_to\_gcf.py
sub\_sample.py
Tensor\_toy.r\_script
text\_to\_bin.py
trim-header.sh
vcf\_2\_map.py

Just the list of scripts for which I need to write up descriptions.
