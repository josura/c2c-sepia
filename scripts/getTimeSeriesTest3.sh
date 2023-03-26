# Loop through all files with the .tsv extension
for file in "*.tsv"; do
  # Extract the prefix name and suffix number
  prefix=$(echo $file | cut -d '-' -f 1)
  suffix=$(echo $file | cut -d '-' -f 2 | cut -d '.' -f 1)
  
  # Extract the desired information and write it to a new file
  awk -v suffix="$suffix" '{if(NR==2) printf("%s\t%s\n", suffix, $5)}' "$file" >> "${prefix}.tsv"
done


