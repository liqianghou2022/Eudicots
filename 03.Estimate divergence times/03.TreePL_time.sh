# Generate IQ-TREE commands for all sampled FASTA files to estimate branch lengths
for f in ../01.random-100/random-100/*.fasta; do  
    echo "/groups/lzu_public/home/houlq2025/software/iqtree-3.0.1-Linux-intel/bin/iqtree3 -s $f -st DNA -m GTR+G -te 01.144-species-4.nwk-qmc.tre -o 143 -T AUTO"
done > run.iqtree.sh

sh run.iqtree.sh

#!/bin/bash

# Iterate over all random_*.tre files in the current directory and generate TreePL configuration files
for treefile in random_*.tre; do
    # Extract the numeric ID, e.g. random_001.tree -> 001
    num=$(echo $treefile | sed -E 's/random_([0-9]+)\.tree/\1/')

    # Name of the configuration file
    cfgfile="${num}.cfg"

    # Write the TreePL configuration file
    cat > "$cfgfile" << EOF
treefile = $treefile
numsites = 1000000
outfile = sample_${num}.out

nthreads = 4
thorough = true
opt = 2
optad = 2
optcvad = 2

mrca = node1 143 124
min = node1 179.9
max = node1 205.0

mrca = node4 41 33
min = node4 102.0
max = node4 112.5

mrca = fossil105 105
max = fossil105 89.8

mrca = fossil18 18
max = fossil18 89.8

mrca = fossil39 39
max = fossil39 93.9

mrca = fossil69 69
max = fossil69 66


cv = true
cvstart = 0.0001
cvstop = 1000
cvmultstep = 0.5
randomcv = true
cviter = 3
cvoutfile = sample_${num}_cv_results.txt

analysis = PL
seed = 12345
EOF

done
