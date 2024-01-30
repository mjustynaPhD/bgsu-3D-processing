#/bin/bash

# for each file in the directory renumbered_pdbs
for file in /tmp/mj/ren/*.zip
do
    filename=$(basename "$file" .zip)
    # unzip, if already unzipped. Save to given name
    unzip -o $file -d /tmp/mj/unzip/$filename
done