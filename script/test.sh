#!/bin/bash

# Organize paths.

# Iterate on cis chromosomes.
for cis in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    echo "***** cis chromosome is $cis!!!!!"
    # Iterate on chromosomes.
    for chromosome in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
    do
        # Determine whether chromosome is cis or trans.
        if [ $chromosome != $cis ]
        then
            # Chromosome is trans.
            echo "chromosome $chromosome is not cis"
            #path="$path_trans_lists/cis"
            #echo path >> trans.txt
        fi
    done
done
