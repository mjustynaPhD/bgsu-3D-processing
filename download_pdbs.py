import argparse
import os
import json
from tqdm import tqdm

SOURCE_RNASOLO = 'https://rnasolo.cs.put.poznan.pl/api/query/structure?'
SOURCE_RCSB = 'https://files.rcsb.org/download/'
SOURCE = SOURCE_RNASOLO
# download pdbs from rcsb.org using pdbid, model and chain
# example: https://files.rcsb.org/download/6z9f_1.pdb


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, default='pdbs', help='Path to output directory')
    parser.add_argument('--json', type=str, default='jsons/il_3.78.json', help='Path to PDB list in JSON format')
    args = parser.parse_args()
    return args

def download(pdbs_chains, pdbs_models, output_dir, skip_existing=True):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for pdb_id, model in tqdm(pdbs_models.items()):
        chains = pdbs_chains[pdb_id]
        for chain in sorted(chains):
            if skip_existing and os.path.exists(f'{output_dir}/{pdb_id}_{model}.pdb'):
                continue
            url = f'{SOURCE}pdbid={pdb_id}&format=PDB&chains={chain}&models={model}'

            # url = f'{SOURCE}{pdb_id}.pdb{model}'
            os.system(f"curl '{url}' --output {output_dir}/{pdb_id}_{model}_{chain}.pdb")
            
            # os.system(f"curl '{url}' --output {output_dir}/{pdb_id}_{model}.pdb")
            # if downloaded file contains "Not Found" string, remove it and download again with cif format
            # with open(f'{output_dir}/{pdb_id}_{model}.pdb') as f:
            #     file_content = f.read()
            #     if 'Not Found' in file_content:
            #         os.remove(f'{output_dir}/{pdb_id}_{model}.pdb')
            #         url = f'{SOURCE}{pdb_id}.cif'
            #         os.system(f"curl '{url}' --output {output_dir}/{pdb_id}.cif")
            #         print(f'Not found: {pdb_id}_{model}')
            #         continue

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
            
    print(f'All pdbs ids in json: {len(pdbs_models)}')
    download(pdb_ids, pdbs_models, args.output_dir)



if __name__ == '__main__':
    main()