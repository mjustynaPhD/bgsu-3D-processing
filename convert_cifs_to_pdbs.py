import os
import argparse
from Bio.PDB import PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from tqdm import tqdm
# Length >=10 nt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, help="Path to input directory")
    parser.add_argument("--output_dir", type=str, help="Path to output directory")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    cifs = os.listdir(args.input_dir)
    cifs = [cif for cif in cifs if cif.endswith(".cif")]
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    for cif in tqdm(cifs):
        basename = cif.replace(".cif", "")
        reader = MMCIFParser()
        structure = reader.get_structure(basename, f'{args.input_dir}/{cif}')
        seq = []
        # change chain id to A
        for model in structure:
            for chain in model:
                if len(chain.id)>1:
                    print(basename, chain.id)
                    chain.id = chain.id[-1]
                for residue in chain:
                    seq.append(residue.get_resname())
                    for atom in residue:
                        # remove altloc
                        atom.set_altloc(" ")
        io = PDBIO()
        io.set_structure(structure)
        io.save(f'{args.output_dir}/{basename}.pdb')
        seq = "".join(seq)
        with open(f'{args.output_dir}/{basename}.seq', 'w') as f:
            f.write(seq)

if __name__=="__main__":
    main()