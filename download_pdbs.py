import argparse
import os
import json
from tqdm import tqdm

SOURCE_RNASOLO = 'https://rnasolo.cs.put.poznan.pl/api/query/structure?'
SOURCE_BGSU = 'http://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/'
SOURCE_RCSB = 'https://files.rcsb.org/download/'
SOURCE = SOURCE_RNASOLO
# download pdbs from rcsb.org using pdbid, model and chain
# example: https://files.rcsb.org/download/6z9f_1.pdb


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, default='pdbs', help='Path to output directory')
    parser.add_argument('--json', type=str, default='jsons/il_3.78.json', help='Path to PDB list in JSON format')
    parser.add_argument('--skip_existing', action='store_true', help='Skip existing files')
    parser.add_argument('--csv_only', action='store_true', help='Download only csv file with pdbs')
    args = parser.parse_args()
    return args

def download_pdbs(pdbs_chains, pdbs_models, output_dir, skip_existing=True):
    print("Downloading pdbs")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for pdb_id, model in tqdm(pdbs_models.items()):
        if pdb_id not in pdbs_chains:
            print(f'Not found: {pdb_id}')
            continue
        chains = pdbs_chains[pdb_id]
        for chain in sorted(chains):
            chs = chain.replace('-', ',')
            if skip_existing and os.path.exists(f'{output_dir}/{pdb_id}_{model}_{chain}.pdb'):
                continue
            url = f'{SOURCE}pdbid={pdb_id}&format=PDB&chains={chs}&models={model}'
            os.system(f"curl '{url}' --output {output_dir}/{pdb_id}_{model}_{chain}.pdb")

def download_csv(pdbs_chains, pdbs_models, output_dir, skip_existing=True):
    print("Downloading csv files")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for pdb_id in tqdm(pdbs_models.keys()):
        if skip_existing and os.path.exists(f'{output_dir}/{pdb_id}.csv'):
            continue
        url = f'{SOURCE_BGSU}{pdb_id}'
        os.system(f"curl '{url}' --output {output_dir}/{pdb_id}.csv")

def load_all_pdb_ids():
    pdb_ids = {}
    with open('data/all_member_pdbids_all__3_300.txt') as f:
        lines = f.readlines()
    
    for l in lines:
        pdb_id = l.split("_")[0]
        chains = l.split("_")[2].strip()
        chains = sorted(chains.split("-"))
        chains = "-".join(chains)
        if pdb_id not in pdb_ids:
            pdb_ids[pdb_id] = []
        pdb_ids[pdb_id].append(chains)
    return pdb_ids

def main():
    args = parse_args()
    # load json file to dict
    with open(args.json) as f:
        json_file = json.load(f)
    pdb_ids = load_all_pdb_ids()
    print(f'All pdbs form txt: {len(pdb_ids)}')
    pdbs_models = {}
    print(f'Internal loops in Json: {len(json_file)}')
    for record in json_file:
        num_instances = record['num_instances']
        num_nucleotides = record['num_nucleotides']
        chainbreak = record['chainbreak']
        alignment = record['alignment']
        for aln in alignment:
            pdbs = record['alignment'][aln]
            pdb_id = pdbs[0].split('|')[0]
            model = pdbs[0].split('|')[1]
            
            pdbs_models[pdb_id] = model
            break  # Get centroids only, i.e. first occurance in json file (class representant)
            
    print(f'All pdbs ids in json: {len(pdbs_models)}')
    if args.csv_only:
        download_csv(pdb_ids, pdbs_models, args.output_dir, skip_existing=args.skip_existing)
    else:
        download_pdbs(pdb_ids, pdbs_models, args.output_dir, skip_existing=args.skip_existing)



if __name__ == '__main__':
    main()