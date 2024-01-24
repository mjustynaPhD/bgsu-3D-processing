# bgsu-3D-processing
Repository to store scripts for structure processing from BGSU RNA 3D motif atlas.


### 1. Downloading the data from RNA 3D motif atlas
Open `http://rna.bgsu.edu/rna3dhub/`. Then go to `RNA 3D Motif Atlas > Home > Internal Loops > Download > json`. The version that I used was [3.78](http://rna.bgsu.edu/rna3dhub/motifs/release/il/3.78)

### 2. Download the structures from RNA Solo

Run `python download_pdbs.py --json data/il_3.78.json --output_dir pdbs`

### 2. Get list of all ILs, HLs with all nucleotides included
For each PDB ID you can download the list of internal loops with ALL nucleotides.

`Note: The json files contains residue numbers that create a loop but do not include the residues that bulge out of the loop. To build the full atomic model we need all residues`
The example link for 6UFG lokks like this: http://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/6UFG

Run `python download_pdbs.py --json data/il_3.78.json --csv_only --output_dir pdbs_csvs`

### 3. Run RNAComposer to generate the structures

#### Create a *.dot file and add GC pairs to stabilize the structure

#### Add GNRA loop to solve the problem of two strands.