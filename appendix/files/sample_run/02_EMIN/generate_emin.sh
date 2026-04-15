#!/bin/bash

# Array of restraint weights
weights=(100 85 75 65 55 45 30 20 15 10 5 1)

# Base content of the input files
base_content='
#200 steps of minimization with explicit solvent and ions and 100.0 kcal/mol-A restraints on protein NOTE 2000 steps would be more typical. This is a short demo!
 &cntrl
     maxcyc=2000, imin=1, ntmin=1, ncyc=250, cut=9.0, igb=0, ntb=1, ntpr=10, ntr=1,
     restraint_wt=PLACEHOLDER, ,restraintmask='\'':1-277'\'',
 &end
END
'

# Create files with different restraint weights
for i in ${!weights[@]}; do
    weight=${weights[$i]}
    filename=$(printf "emin%d.in" $((i+1)))
    content=$(echo "$base_content" | sed "s/PLACEHOLDER/$weight/")
    echo "$content" > "$filename"
done

