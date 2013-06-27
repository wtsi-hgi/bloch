#!/bin/bash

for i in {1..14}; do
 echo "i=$i"
 python -m cProfile -o ./TreeBuildingProfiles/g$i_tree.profile PhasingFinal2.py -t /Users/mp18/Documents/bloch/Data/genotype_$i.txt -o /Users/mp18/Documents/bloch/Data/SerialisedTrees/g$i
done
