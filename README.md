# bgsu-3D-processing
Repository to store scripts for structure processing from BGSU RNA 3D motif atlas.


### Downloading the data from RNA 3D motif atlas
Open `http://rna.bgsu.edu/rna3dhub/`. Then go to `RNA 3D Motif Atlas > Home > Internal Loops > Download > json`. The version that I used was [3.78](http://rna.bgsu.edu/rna3dhub/motifs/release/il/3.78)

### Download the structures from RNA Solo

Run `python download_pdbs.py --json data/il_3.78.json --output_dir pdbs`

### Get list of all ILs, HLs with all nucleotides included
For each PDB ID you can download the list of internal loops with ALL nucleotides.

`Note: The json files contains residue numbers that create a loop but do not include the residues that bulge out of the loop. To build the full atomic model we need all residues`
The example link for 6UFG lokks like this: http://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/6UFG

Run `python download_pdbs.py --json data/il_3.78.json --csv_only --output_dir pdbs_csvs`

### Run RNAComposer to generate the structures

#### Renumber the motifs that you want to use and generate *.dot files
For this purpose generate the commands to run RNAqua or RNA-tools. The command should look like this:
`process_for_RNAqua.py --directory pdbs_csv/ --output_dir RNAtools --rna-tools`

The parameter `--rna-tools` is used to generate the command for RNA-tools. If you want to use RNAqua, just remove the parameter.
For my usage I used RNA-tools option. However, this required modified version of RNA-tools script. The original one after donwload has a bug, which was fixed by me. Another problem is that the PDB files remain big and RNA-Composer process these files significantly longer.

You can use RNAaqua, but for some files I had a problem with the output. In some cases the output zip file was corrupted and could not be used. The were two main reasons for this:
1. The files contained some non-standard characters in the file names. This was fixed by replacing '-' with '_'.
2. Some files names contain incorrect chain id, e.g. 1XYZ_C, while the chain id in file is A.

This command will also generate *.dot files, required for RNA-Composer, with GNRA loop added to the structure and additional GC pairs.

The output processed for RNAaqua is in zip format, so then use `./unzip_all_files.sh` to unzip all files.

#### Generate sbs files
Once the renumbered pdb files are generated you can run script `prepare_sbs.sh` to generate sbs files. The script will generate sbs files for each pdb file in the directory. The *.sbs and *.dot files are required for RNA-Composer in the next step.

#### (Optional) Transfer files to the server
You can use `rsync` to transfer the generated sbs files and dot files to the server. The command should look like this:
`rsync -av sbs_files/ user@server:/path/to/sbs_files`

#### (Optional) Fix file names
Note that RNA-Composer does not accept '-' in the file names. You need to replace '-' with '_'. For this purpose I used the following command:
`rename 's/-/_/g' *` or `python fix_names.py <path>`

#### Run RNA-Composer
You can now run RNA-Composer on server. Run the following command:
`./run_composer.sh`