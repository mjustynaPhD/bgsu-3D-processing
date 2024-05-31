import argparse
import json
import os
from tqdm import tqdm

DOWNLOAD_URL = "http://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord="
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=str, default="jsons/il_3.78.json", help="Path to PDB list in JSON format")
    parser.add_argument("--centroid-only", action="store_true", help="Use centroid only")
    parser.add_argument("--output_dir", type=str, default="pdbs", help="Path to output directory")
    args = parser.parse_args()
    return args

def extract_motif_only(motif_file):
    with open(motif_file) as f:
        lines = f.readlines()
    trim_id = None
    for i, line in enumerate(lines):
        if i>0 and line.startswith("data_view"):
            trim_id = i
            break
    with open(motif_file, "w") as f:
        lines = lines[:trim_id]
        lines = [l.replace("?", "0") for l in lines]
        f.writelines(lines)   


def main():
    args = parse_args()
    with open(args.json) as f:
        json_file = json.load(f)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for motif in tqdm(json_file):
        for align in motif['alignment'].keys():
            url = f'{DOWNLOAD_URL}{align}'
            # download url and save in output_dir
            os.system(f"curl '{url}' --output {args.output_dir}/{align}.cif")
            extract_motif_only(f"{args.output_dir}/{align}.cif")
            if args.centroid_only:
                break
            


    pass

if __name__=="__main__":
    main()