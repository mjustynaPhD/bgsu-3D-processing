import argparse
import os
import json
from tqdm import tqdm

SOURCE_RNASOLO = 'https://rnasolo.cs.put.poznan.pl/api/query/structure?'
SOURCE_RCSB = 'https://files.rcsb.org/download/'
SOURCE = SOURCE_RCSB
# download pdbs from rcsb.org using pdbid, model and chain
# example: https://files.rcsb.org/download/6z9f_1.pdb


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, default='pdbs', help='Path to output directory')
    parser.add_argument('--json', type=str, default='jsons/il_3.78.json', help='Path to PDB list in JSON format')
    args = parser.parse_args()
    return args

def download(pdbs, output_dir, skip_existing=True):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for pdb in tqdm(pdbs):
        
        pdb_id, model, chain = pdb
        if skip_existing and os.path.exists(f'{output_dir}/{pdb_id}_{model}.pdb'):
            continue
        # url = f'{SOURCE}pdbid={pdb_id}&format=PDB&chains={chain}&models={model}'
        url = f'{SOURCE}{pdb_id}.pdb{model}'
        # os.system(f"curl '{url}' --output {output_dir}/{pdb_id}_{model}_{chain}.pdb")
        
        os.system(f"curl '{url}' --output {output_dir}/{pdb_id}_{model}.pdb")
        # if downloaded file contains "Not Found" string, remove it and download again with cif format
        with open(f'{output_dir}/{pdb_id}_{model}.pdb') as f:
            file_content = f.read()
            if 'Not Found' in file_content:
                os.remove(f'{output_dir}/{pdb_id}_{model}.pdb')
                url = f'{SOURCE}{pdb_id}.cif'
                os.system(f"curl '{url}' --output {output_dir}/{pdb_id}.cif")
                print(f'Not found: {pdb_id}_{model}')
                continue

def main():
    args = parse_args()
    # load json file to dict
    with open(args.json) as f:
        json_file = json.load(f)
    
    pdb_to_download = []
    print(len(json_file))
    for record in json_file:
        num_instances = record['num_instances']
        num_nucleotides = record['num_nucleotides']
        chainbreak = record['chainbreak']
        alignment = record['alignment']
        for aln in alignment:
            pdbs = record['alignment'][aln]
            pdb_id = pdbs[0].split('|')[0]
            model = pdbs[0].split('|')[1]
            chain = pdbs[0].split('|')[2]
            pdb = (pdb_id, model, chain)
            pdb_to_download.append(pdb)
    print(f'All instances: {len(pdb_to_download)}')
    pdb_to_download = list(set(pdb_to_download)) # remove duplicates
    print(f'Unique pdbs with chains: {len(pdb_to_download)}')
    download(pdb_to_download, args.output_dir)



if __name__ == '__main__':
    main()