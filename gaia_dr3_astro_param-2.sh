#!/bin/bash

filelist="$PWD/gaiafiles-2.txt"
logfile="$PWD/gaiadr3_astro_param_2_$(date +"%m-%d-%Y").log"

for file in $(cat $filelist);
do
now=$(date)
echo "# $now Start downloading $file." | tee -a "$logfile"
wget http://cdn.gea.esac.esa.int/Gaia/gdr3/Astrophysical_parameters/astrophysical_parameters/"$file" | tee -a "$logfile"

now=$(date)
echo "# $now Download finished, start extracting $file." | tee -a "$logfile"
gzip -dk "$PWD/$file"

now=$(date)
echo "# $now Extraction finished, start normalizing file content (${file:0:41})." | tee -a "$logfile"
awk '!/^#/' "$PWD/${file:0:41}" > "$PWD/tmp_${file:0:41}"
#python3 gaia_list_filter.py "$PWD/tmp_${file:0:28}" | tee -a "$logfile"
python3 ~/Documents/dev/github/gaiadr3_double_stars/gaia_dr3_astro_param.py "$PWD/tmp_${file:0:41}" "$PWD/gaia_doublestar_sourceid.txt"
mv "$PWD/result.csv" "$PWD/final_${file:0:41}"

#now=$(date)
#echo "# $now Normalization finished, start searching for double stars." | tee -a "$logfile"
#python3 search_doubles_gaiadr3_v2.py "$PWD/final_${file:0:28}" > "$PWD/doubles_${file:0:28}" | tee -a "$logfile"

now=$(date)
echo "# $now Identification of double stars finished, results can be found in file doubles_"${file:0:41}". Start removing temporary files." | tee -a "$logfile"
rm "$PWD/$file" "$PWD/${file:0:41}" "$PWD/tmp_${file:0:41}"

now=$(date)
echo "# $now Temporary files removed." | tee -a "$logfile"
done