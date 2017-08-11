tempDir=$PWD/temperate/merged/
tropDir=$PWD/tropical/merged/

## Get names of bam files for each population
find $tempDir -name *.bam > temperateList.txt
find $tropDir -name *.bam > tropicalList.txt

## Make inbreeding coefficient files
## Assume all lines fully inbred
for i in $(seq 30)
do
    echo 1.0	
done > temperate.indF

for i in $(seq 30)
do
    echo 1.0
done > tropical.indF

angsd=$PWD/angsd/angsd
reference=$PWD/Zea_mays.AGPv3.31.dna.genome.fa
minQ=20
minMapQ=30

nInd=$( wc -l temperateList.txt | cut -f1 -d" " )
minProportion=0.7
minInd=$( printf "%.0f" $(echo "scale=2;$minProportion*$nInd" | bc))

mkdir -p saf_windows_split/
rm -f saf_windows_split/*

#split -l 600 allFst_windows.txt $PWD/saf_windows_split/
for chr in $( seq 10 )
do
    grep -P "^$chr:" allFst_windows.txt > saf_windows_split/$chr
done

for pop in temperate tropical
do
	mkdir -p $PWD/${pop}_saf/
	rm -f $PWD/${pop}_saf/*
	
	angsdCommand="nohup $angsd \
          -bam ${pop}List.txt \
	  -minQ $minQ \
	  -minMapQ $minMapQ \
	  -baq 1 \
	  -ref $reference \
	  -anc $reference \
	  -doCounts 1 \
	  -doQsDist 1 \
	  -doDepth 1 \
          -rf {} \
          -P 4 \
	  -GL 1 \
	  -doSaf 2 \
	  -doMaf 1 \
	  -doMajorMinor 1 \
	  -indF ${pop}.indF \
	  -nInd $nInd \
	  -minInd $minInd \
	  -out $PWD/${pop}_saf/${pop}_saf.{/} &> $PWD/${pop}_saf/${pop}_saf.{/}.log &"

	ls saf_windows_split/* | parallel $angsdCommand
done


