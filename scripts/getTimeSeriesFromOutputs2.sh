#!/usr/bin/env bash



# set the directory containing the input files
input_dir=../outputs

# set the output directory where the output files will be saved
output_dir=../outputsTimeSeries


# loop over each unique prefix name in the input directory
for prefix in $(ls $input_dir | cut -d '-' -f 1 | sort -u); do

  # initialize the output file with the feature column header
  echo -e "Feature\t$prefix" > $output_dir/${prefix}_output.tsv

 # loop over each file with the current prefix
  for file in $(ls $input_dir | grep "^$prefix"); do
      out_file="$output_dir/${file%.*}.out"
      awk -F'\t' 'FNR==1 {next} {print substr(FILENAME, index(FILENAME, "--")+2), $5}' $file > $out_file
  done
done