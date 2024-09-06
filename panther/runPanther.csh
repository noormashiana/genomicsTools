#!/bin/tcsh

# Incomplete, paths are still hardcoded for testing

set species = $1
set num_threads = 8

set start_time = `date +%s`

set path_to_fastas = /pscratch/sd/m/mashiana/Vertebrates/fastas
set output_directory = /pscratch/sd/m/mashiana/Vertebrates/panther

mkdir -p ${output_directory}/tables
mkdir -p ${output_directory}/error

cd /pscratch/sd/m/mashiana/panther/pantherScore2.2
source panther.cshrc

mkdir /tmp/${species}

./pantherScore2.2.pl -l /pscratch/sd/m/mashiana/panther/target/famlib/rel/PANTHER18.0_altVersion/hmmscoring/PANTHER18.0 -D B -V -i ${path_to_fastas}/${species}.fasta -o ${output_directory}/tables/${species}.tsv -e ${output_directory}/error/${species}.error -n -c ${num_threads} -T /tmp/${species} 


if ($status == 0) then
    set end_time = `date +%s`

    @ elapsed_time = $end_time - $start_time

    echo "PANTHER scoring for ${species} completed successfully." >> ${output_directory}/error/${species}.log
    echo "Elapsed time: ${elapsed_time} seconds." >> ${output_directory}/error/${species}.log
    
    # Remove the temporary directory
    rm -rf /tmp/${species}
else
    echo "PANTHER scoring for ${species} failed." >> -a ${output_directory}/error/${species}.log
endif
