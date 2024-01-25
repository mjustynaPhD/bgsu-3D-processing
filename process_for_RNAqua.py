import os
import argparse
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description='Process RNAqua output for further analysis')
    parser.add_argument('-i', '--input', type=str, help='Input file name')
    parser.add_argument('-d', '--directory', type=str, help='Process all files in the directory')
    parser.add_argument('-o', '--output_dir', type=str, help='Output directory')
    # parser.add_argument('-il', '--include_il', action='store_true', help='Include IL in the output')
    return parser.parse_args()

def process_file(input_file, output_dir, il:bool=True, hl:bool=False, j3:bool=False):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if il:
        lines = [l for l in lines if 'IL' in l]
    commands = []
    for l in lines:
        # l = lines[212]
        l = l.strip().split('","')
        name = l[0].replace('"', '')
        residues = l[1].split(',')
        begin_end = l[2].replace('"', '').split(',')
        begin_end = [int(x) for x in begin_end]
        pdb_name = "_".join(residues[0].split('|')[:3])
        pdb_name2 = "_".join(residues[-2].split('|')[:3])
        if pdb_name != pdb_name2:
            second_chain = pdb_name2[-1]
            pdb_name = f'{pdb_name}-{second_chain}'
        seqs, rev_seqs, res_ids, rev_res_ids = extract_sequence(residues, begin_end)  # TODO: save the alignments with the sequences in fasta format
        cmd = get_command(pdb_name, res_ids)
        cmd2 = get_command(pdb_name, rev_res_ids, rev=True)
        commands.append(cmd)
        commands.append(cmd2)
    return commands

def process_directory(input_dir, output_dir):
    dir_files = os.listdir(input_dir)
    commands = []
    for file in dir_files:
        if file.endswith('.csv'):
            cmds = process_file(os.path.join(input_dir, file), output_dir)
            commands.extend(cmds)
    return commands

def extract_sequence(residues, begin_end):
    sequence = [r.split('|')[3] for r in residues]
    res_nums = [r.split('|')[4] for r in residues]
    chains = [r.split('|')[2] for r in residues]
    seq = ''.join(sequence)
    seqs, bgend_pairs = split_sequence(seq, begin_end)
    res_ids, rev_res_ids = get_residue_ids(res_nums, chains, bgend_pairs)
    rev_seqs = [s for s in seqs[::-1]]
    return seqs, rev_seqs, res_ids, rev_res_ids

def split_sequence(seq, begin_end):
    # get indeces where begin_end is 1
    begin_end = np.array(begin_end)
    begin_end = np.where(begin_end == 1)[0]
    begin_end = begin_end.reshape(-1, 2)
    seqs = []
    for be in begin_end:
        seqs.append(seq[be[0]:be[1]+1])
    return seqs, begin_end

def get_residue_ids(res_nums, chains, bgend_pairs):
    strands = []
    for a, b in bgend_pairs:
        ch = chains[a]
        if len(ch) == 2:
            ch = ch[1]
        strands.append(f"{ch}_{res_nums[a]},{(b-a)+1}")
    joined = "|".join(strands)
    rev_joined = "|".join(strands[::-1])
    return joined, rev_joined

def get_command(file_name, ids, rev:bool=False):
    if rev:
        command = f'docker-compose run --rm --entrypoint ./rnaqua rnaqua --command RENUMERATED-3D --single-model-file-path /tmp/mj/{file_name}.pdb --alignment "{file_name}.pdb:{ids}" --output-file-path /tmp/mj/{file_name}_m_rev.zip'
    else:
        command = f'docker-compose run --rm --entrypoint ./rnaqua rnaqua --command RENUMERATED-3D --single-model-file-path /tmp/mj/{file_name}.pdb --alignment "{file_name}.pdb:{ids}" --output-file-path /tmp/mj/{file_name}_m.zip'
    return command

def write_commands(commands, output_dir):
    with open(output_dir, 'w') as f:
        for c in commands:
            f.write(f"{c}\n")

def main():
    args = parse_args()
    if args.input is not None and args.directory is not None:
        print('Please specify either input file or directory')
        return
    
    if args.input is not None:
        cmds = process_file(args.input, args.output_dir)
    elif args.directory is not None:
        cmds = process_directory(args.directory, args.output_dir)
    else:
        print('Please specify either input file or directory')
        return
    assert args.output_dir is not None
    assert len(cmds) > 0
    print(f'Number of commands: {len(cmds)}')
    write_commands(cmds, args.output_dir)

if __name__ == '__main__':
    main()