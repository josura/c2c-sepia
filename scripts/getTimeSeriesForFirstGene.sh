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
    echo "getting first gene value for file $file"

    # get the feature name from the file name
    feature=$(echo $file | cut -d '-' -f 3 | cut -d '.' -f 1)

    # get the fifth value from the current file
    value=$(awk 'NR==5 {print $0}' $input_dir/$file)

    echo -n $file | cut -d '-' -f 3 | cut -d '.' -f 1 | tr -d $'\n' >> $output_dir/${prefix}_output.tsv

    # append the feature and value to the output file
    cat $input_dir/$file | awk '{if(NR==2) printf("\t%s\n",$5)}'  >> $output_dir/${prefix}_output.tsv
    #echo -e "$feature\t$value" >> $output_dir/${prefix}_output.tsv

  done

done
