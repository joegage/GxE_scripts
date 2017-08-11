# ## Get names of bam files for each population
# find $tempDir -name *.bam > temperateList.txt
# find $tropDir -name *.bam > tropicalList.txt
# 
# ## Make inbreeding coefficient files
# ## Assume all lines fully inbred
# for i in $(seq 30)
# do
#     echo 1.0
# done > temperate.indF
# 
# for i in $(seq 30)
# do
#     echo 1.0
# done > tropical.indF
# 
# cat temperateList.txt tropicalList.txt > allList.txt
# cat temperate.indF tropical.indF > all.indF

angsd=$PWD/angsd/angsd
reference=$PWD/Zea_mays.AGPv3.31.dna.genome.fa
minQ=20
minMapQ=30

nInd=30
minProportion=0.7
minInd=$( printf "%.0f" $(echo "scale=2;$minProportion*$nInd" | bc))

mkdir -p nd_windows_split/
rm -f nd_windows_split/*
split -l 600 allFst_windows.txt $PWD/nd_windows_split/
#for chr in $(seq 10)
#do
# 	grep "${chr}:" allFst_windows.txt > $PWD/nd_windows_split/chr${chr}
#done

# for pop in temperate tropical
# do
# 	for chr in $(seq 10)
# 	do
# 		mkdir -p $PWD/${pop}_nd/
		#angsdCommand="nohup $angsd -bam ${pop}List.txt -minQ $minQ -minMapQ $minMapQ -baq 1 -ref $reference -anc $reference -doCounts 1 -doQsDist 1 -doDepth 1 -rf {} -P 32 -GL 1 -doSaf 2 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-2 -doThetas 1 -pest ${pop}_saf_wg/${pop}_saf.sfs -indF ${pop}.indF -nInd $nInd -minInd $minInd -out $PWD/${pop}_nd/${pop}_nd.{/} &> $PWD/${pop}_nd/${pop}_nd.{/}.log &"
		

for pop in temperate tropical
do
	angsdCommand="nohup $angsd -bam ${pop}List.txt -minQ $minQ -minMapQ $minMapQ -baq 1 -ref $reference -anc $reference -doCounts 1 -doQsDist 1 -doDepth 1 -rf {} -GL 1 -doSaf 2 -doMaf 1 -doMajorMinor 1 -doThetas 1 -pest ${pop}_saf/${pop}_saf.sfs -indF ${pop}.indF -nInd $nInd -minInd $minInd -out $PWD/${pop}_nd/${pop}_nd.{/} &> $PWD/${pop}_nd/${pop}_nd.{/}.log &"
	
	mkdir -p ${pop}_nd
	rm -f ${pop}_nd/*
	ls nd_windows_split/* | parallel $angsdCommand
done
