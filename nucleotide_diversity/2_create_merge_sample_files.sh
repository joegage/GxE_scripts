## Some of the 30 temperate and 30 tropical lines have multiple
## samples and therefore multiple bam files.  This will create
## text files with a list of bam files for each genotype.

mkdir -p ./temperate/sample_lists/
while read geno
do
    find $PWD/temperate/ -maxdepth 1 -name $geno*.bam > ./temperate/sample_lists/${geno}_samples.txt
done < temperate_inbreds.txt

mkdir -p ./tropical/sample_lists/
while read geno
do
    find $PWD/tropical/ -maxdepth 1 -name $geno*.bam > ./tropical/sample_lists/${geno}_samples.txt
done < tropical_inbreds.txt

