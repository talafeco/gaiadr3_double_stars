#!/bin/bash

filelist="$PWD/gaiafiles.txt"
logfile="$PWD/gaiadr3_doubles_$(date +"%m-%d-%Y").log"

for file in $(cat $filelist);
do
now=$(date)
echo "# $now Start downloading $file." | tee -a "$logfile"
wget https://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/"$file" | tee -a "$logfile"

now=$(date)
echo "# $now Download finished, start extracting $file." | tee -a "$logfile"
gzip -dk "$PWD/$file"

now=$(date)
echo "# $now Extraction finished, start normalizing file content (${file:0:28})." | tee -a "$logfile"
awk '!/^#/' "$PWD/${file:0:28}" > "$PWD/tmp_${file:0:28}"
python3 ~/Git/gaiadr3_double_stars/tools/gaia_list_filter_v2.py "$PWD/tmp_${file:0:28}" | tee -a "$logfile"
mv "$PWD/result.csv" "$PWD/final_${file:0:28}"

#now=$(date)
#echo "# $now Normalization finished, start searching for double stars." | tee -a "$logfile"
#python3 search_doubles_gaiadr3_v2.py "$PWD/final_${file:0:28}" > "$PWD/doubles_${file:0:28}" | tee -a "$logfile"

now=$(date)
echo "# $now Identification of double stars finished, results can be found in file doubles_"${file:0:28}". Start removing temporary files." | tee -a "$logfile"
rm "$PWD/$file" "$PWD/${file:0:28}" "$PWD/tmp_${file:0:28}"

now=$(date)
echo "# $now Temporary files removed." | tee -a "$logfile"
done