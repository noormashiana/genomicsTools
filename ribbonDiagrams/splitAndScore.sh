#!/bin/bash

# better for scoring fewer bigger fastas
# requires hmmscan, seqkit, seqkit I used modules one

start_time=$(date +%s)

num_split="$1"
source_file="$2"
hmmlib_path="$3"
split_directory="$4"

destination_directory="$5"

confusion_directory="$6"

process_file() {
    local file="$1"
    local species_name=$(basename "$file")
    local output_file="${destination_directory}/${species_name}"
    echo "$species_name"
    echo "$output_file"
    echo "$hmmlib_path"
    echo "$file"
    /global/homes/m/mashiana/tools/hmmer3.1b2/bin/hmmscan --cpu 1 --tblout "$output_file" "$hmmlib_path" "$file" &> /dev/null
    awk -v species_name="$species_name" '!/^#/ && !seen[$3]++ {print $3"\t"species_name"\t"$1"\t"$5"\t"$6}' "$output_file" > "$confusion_directory/$species_name"
}

export -f process_file
export destination_directory
export hmmlib_path
export confusion_directory

mkdir -p "$destination_directory"
mkdir -p "$confusion_directory"
mkdir -p "$split_directory"

# change this to seqkit
/global/cfs/cdirs/plantbox/tools/perlmutter/modules/seqkit/2.3.1/bin/seqkit split -p "$num_split" -O "$split_directory" "$source_file"

find "$split_directory" -type f -name '*' | parallel --gnu -j "$num_split" process_file {} {/}

end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Total execution time: $execution_time seconds"