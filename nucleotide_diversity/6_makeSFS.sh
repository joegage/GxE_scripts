for pop in temperate tropical
do

    echo "Catting all the saf files..."
    ./angsd/misc/realSFS cat ${pop}_saf/${pop}_saf.1.saf.idx ${pop}_saf/${pop}_saf.2.saf.idx ${pop}_saf/${pop}_saf.3.saf.idx ${pop}_saf/${pop}_saf.4.saf.idx ${pop}_saf/${pop}_saf.5.saf.idx ${pop}_saf/${pop}_saf.6.saf.idx ${pop}_saf/${pop}_saf.7.saf.idx ${pop}_saf/${pop}_saf.8.saf.idx ${pop}_saf/${pop}_saf.9.saf.idx ${pop}_saf/${pop}_saf.10.saf.idx -outnames ${pop}_saf/${pop}_saf
    
    echo "Will the real SFS please stand up?"
    ./angsd/misc/realSFS ${pop}_saf/${pop}_saf.saf.idx -P 12 > ${pop}_saf/${pop}_saf.sfs
done
