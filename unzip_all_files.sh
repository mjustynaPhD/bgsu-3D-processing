#/bin/bash

# for each file in the directory renumbered_pdbs
for file in /tmp/mj/ren/*.zip
do
    unzip $file -d /tmp/mj/ren
done