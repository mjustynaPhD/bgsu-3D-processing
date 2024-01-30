#!/bin/bash

# # for each file in the directory renumbered_pdbs
# for file in unzip/*.zip
# do
#     # get the filename without the extension
#     filename=$(basename "$file" .zip)
#     # run prepare_sbs.py on the file
#     # python prepare_sbs.py $file
#     # move the output file to the directory sbs_pdbs
#     # mv $filename.sbs.pdb sbs_pdbs
#     echo $filename
#     echo -e "1;1\nL1" > sbs_files/$filename.sbs; cat unzip/$filename/*.pdb >> sbs_files/$filename.sbs
# done

# for each file in the directory renumbered_pdbs
for file in renumbered_pdbs/*.pdb
do
    # get the filename without the extension
    filename=$(basename "$file" .pdb)
    echo $filename
    echo -e "1;1\nL1" > sbs_files/$filename.sbs; cat renumbered_pdbs/$filename.pdb >> sbs_files/$filename.sbs
done