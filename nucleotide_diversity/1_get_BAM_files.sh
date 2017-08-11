## Get file list of the hapmap directory
if [ ! -e fileList.txt ]
  then
    ils -r --bundle /iplant/home/shared/panzea/hapmap3/bam/ > fileList.txt
fi

## Find the necessary files
grep -f temperate_inbreds.txt fileList.txt > temperate_file_paths.txt
grep -f tropical_inbreds.txt fileList.txt > tropical_file_paths.txt

## Trim annoying bits off the lines with file paths
sed -i 's/Bundle file: //' temperate_file_paths.txt
sed -i 's/Bundle file: //' tropical_file_paths.txt

mkdir -p temperate
mkdir -p tropical

echo -e "Downloading files: "

## Download temperate files
while read file
do
    if [ -e ./temperate/$file ]
      then
        echo "$file already downloaded"
      else
        echo "Downloading $file"
        iget $file ./temperate/
    fi
done < temperate_file_paths.txt

## Remove one of the B73 samples - no index file, seems to be corrupted or
## broken in some way...
rm ./temperate/B73_C11VBACXX_ACAGTG.bam

## Download tropical files
while read file
do
    if [ -e ./tropical/$file ]
      then
        echo "$file already downloaded"
      else
        echo "Downloading $file"
        iget $file ./tropical/
    fi
done < tropical_file_paths.txt
