#!/bin/bash

folder_path=$1

max_lines=0

for file in "$folder_path"/*.tsv; do
    if [[ -f "$file" ]]; then
        lines=$(wc -l < "$file")
        if (( lines > max_lines )); then
            max_lines=$lines
        fi
    fi
done

echo "Maximum number of lines: $max_lines"
