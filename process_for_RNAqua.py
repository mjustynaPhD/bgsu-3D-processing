import os
import argparse
import json
import numpy as np
except_pdbs = {'5TBW_1_1':'5TBW_1_1-4',
               '5TBW_1_4':'5TBW_1_1-4',
               '7JRS_1_B':'7JRS_1_A-B',
               '1U6B_1_B-C': '1U6B_1_B',
               '3P59_1_E-F': '3P59_1_A-B-C-D-E-F-G-H'}

def parse_args():
    parser = argparse.ArgumentParser(description='Process RNAqua output for further analysis')
    parser.add_argument('-i', '--input', type=str, help='Input file name')
    parser.add_argument('-d', '--directory', type=str, help='Process all files in the directory')
    parser.add_argument('-o', '--output_dir', type=str, help='Output directory')
    parser.add_argument('--rna-tools', action='store_true', help='Prepare commands for rna-tools package')
    # parser.add_argument('-il', '--include_il', action='store_true', help='Include IL in the output')
    return parser.parse_args()

def process_file(input_file, output_dir, il:bool=True, hl:bool=False, j3:bool=False, centroids:list=None, rna_tools:bool=False):
    """
    This procedure iterates over a csv file and prepares commands for each instance in the csv file.
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if il:  # include internal loops
        lines = [l for l in lines if 'IL' in l]
    if centroids is not None:  # if centroid is specified, include only this centroid. Otherwise, include all
        lines = [l for l in lines if l.split(',')[0].replace('"', '') in centroids]

    commands = []
    for l in lines:
        # l = lines[212]
        l = l.strip().split('","')
        name = l[0].replace('"', '')
        residues = l[1].split(',')
        begin_end = l[2].replace('"', '').split(',')
        begin_end = [int(x) for x in begin_end]  # extract chainbraker positions
        pdb_name = "_".join(residues[0].split('|')[:3])
        pdb_name2 = "_".join(residues[-2].split('|')[:3])
        if pdb_name != pdb_name2:
            second_chain = pdb_name2[-1]
            pdb_name = f'{pdb_name}-{second_chain}'
        if pdb_name in except_pdbs:
            pdb_name = except_pdbs[pdb_name]
        pdb_name = pdb_name.replace('-', '_') # replace '-' with '_' in multi-chains, otherwise RNA-Composer will not work
        seqs, rev_seqs, res_ids, rev_res_ids = extract_sequence(
                                                        residues,
                                                        begin_end,
                                                        rna_tools=rna_tools)
        cmd = get_command(pdb_name, res_ids, loop_name = name, rna_tools=rna_tools)
        cmd2 = get_command(pdb_name, rev_res_ids, loop_name = name, rev=True, rna_tools=rna_tools)
        write_dot_bracket(f"{pdb_name}_{name}", seqs, begin_end, output_dir)
        write_dot_bracket(f"{pdb_name}_{name}", rev_seqs, begin_end, output_dir, rev=True)

        commands.append(cmd)
        commands.append(cmd2)
    return commands

def process_directory(input_dir, output_dir, centroids:list=None, rna_tools:bool=False):
    """
    Iterate over all csv files for each pdb and prepare commands for each instance in every csv file.
    """
    dir_files = os.listdir(input_dir)
    if centroids is not None:
        cetnroids_pdbs = [f.split('_')[1] for f in centroids]
    else:
        cetnroids_pdbs = None
    commands = []
    for file in dir_files:
        pdb_name = file.replace(".csv", "")
        if file.endswith('.csv') and pdb_name in cetnroids_pdbs:
            cmds = process_file(os.path.join(input_dir, file), output_dir, centroids=centroids, rna_tools=rna_tools)
            commands.extend(cmds)
    return commands

def extract_sequence(residues, begin_end, rna_tools:bool=False):
    """
    This procedure prepares sequences and residue ids for each instance in the csv file.
    The sequences are returned as a list of strings, where each string is a sequence for one instance, e.g.
    seqs = ['CACGGCG', 'CGG']
    The residue ids are returned as a string, where each instance is separated by a comma. Each instance is
    prepared to be processed by tools such as RNAqua or rna-tools. The example below shows the format for
    rna-tools:
    res_ids = 'B:60-66>B:1-7,B:78-80>B:1-3'
    rev_res_ids = 'B:78-80>B:1-3,B:60-66>B:1-7'

    Args:
        residues (list): list of residues from the csv file, e.g. ['5U3G|1|B|C|60', '5U3G|1|B|C|61', ...]
        begin_end (list): list of chainbraker positions
        rna_tools (bool): if True, prepare residue ids for rna-tools, otherwise prepare for RNAqua

    """
    sequence = [r.split('|')[3] for r in residues]
    res_nums = [r.split('|')[4] for r in residues]
    chains = [r.split('|')[2] for r in residues]
    seq = ''.join(sequence)
    seqs, bgend_pairs = split_sequence(seq, begin_end)
    if rna_tools:
        res_ids, rev_res_ids = get_residue_ids_rna_tools(
                                                        res_nums,
                                                        chains,
                                                        bgend_pairs)
    else:
        res_ids, rev_res_ids = get_residue_ids_aqua(res_nums, chains, bgend_pairs)
    rev_seqs = [s for s in seqs[::-1]]
    return seqs, rev_seqs, res_ids, rev_res_ids

def split_sequence(seq, begin_end):
    """
    Split sequence into subsequences based on chainbraker positions.
    The indeces, where begin_end is 1, are used to split the nucleotide sequence.
    Args:
        seq (str): nucleotide sequence
        begin_end (list): list of chainbraker positions
    Returns:
        seqs (list): list of subsequences
        begin_end (list): list of begin and end positions for each subsequence
    """
    begin_end = np.array(begin_end)
    begin_end = np.where(begin_end == 1)[0]
    begin_end = begin_end.reshape(-1, 2)
    seqs = []
    for be in begin_end:
        seqs.append(seq[be[0]:be[1]+1])
    return seqs, begin_end

def get_residue_ids_rna_tools(res_nums, chains, bgend_pairs):
    """
    Prepare residue ids in the format for rna-tools package.
    """
    strands = []
    prev_end = 1
    for a, b in bgend_pairs:
        diff = int(res_nums[b]) - int(res_nums[a])
        ch = chains[a]
        if (len(ch) == 2 and 'A' in ch) or ch.isdigit():
            ch = 'A'
        elif len(ch) == 2:
            ch = ch[1]
        strands.append(f"{ch}:{res_nums[a]}-{res_nums[b]}>{ch}:{prev_end}-{prev_end+diff}")  # A:201-208>A:1-8
        if prev_end == 1:
            prev_end += diff + 1

    joined = ",".join(strands)
    rev_joined = ",".join(strands[::-1])
    return joined, rev_joined

def get_residue_ids_aqua(res_nums, chains, bgend_pairs):
    """
    Prepare residue ids in the format for RNAqua package.
    """
    strands = []
    for a, b in bgend_pairs:
        ch = chains[a]
        if len(ch) == 2:
            ch = ch[1]
        strands.append(f"{ch}_{res_nums[a]},{(b-a)+1}") # 'A_24,7'
    joined = "|".join(strands)
    rev_joined = "|".join(strands[::-1])
    return joined, rev_joined

def get_command(file_name, ids, loop_name, rev:bool=False, rna_tools:bool=False):
    if rna_tools:
        command = get_rna_tools_command(file_name, ids, loop_name=loop_name, rev=rev)
    else:
        command = get_aqua_command(file_name, ids, loop_name=loop_name, rev=rev)
    return command

def get_aqua_command(file_name, ids, loop_name:str, rev:bool=False):
    if rev:
        command = f'docker-compose run --rm --entrypoint ./rnaqua rnaqua --command RENUMERATED-3D --single-model-file-path /tmp/mj/pdbs/{file_name}.pdb --alignment "{file_name}.pdb:{ids}" --output-file-path /tmp/mj/ren/{file_name}_{loop_name}_rev.zip'
    else:
        command = f'docker-compose run --rm --entrypoint ./rnaqua rnaqua --command RENUMERATED-3D --single-model-file-path /tmp/mj/pdbs/{file_name}.pdb --alignment "{file_name}.pdb:{ids}" --output-file-path /tmp/mj/ren/{file_name}_{loop_name}.zip'
    return command

def get_rna_tools_command(file_name, ids, loop_name:str, rev:bool=False):
    """
    Prepare commands as in example below:
    `rna_pdb_tools.py --edit 'A:201-208>A:1-8,A:188-195>A:9-16' <pdb_file>.pdb >pdb_file_edited_<num>.pdb`

    """
    if rev:
        command = f"rna_pdb_tools.py --edit '{ids}' pdbs/{file_name}.pdb > renumbered_pdbs/{file_name}_{loop_name}_rev.pdb"
    else:
        command = f"rna_pdb_tools.py --edit '{ids}' pdbs/{file_name}.pdb > renumbered_pdbs/{file_name}_{loop_name}.pdb"
    return command

def write_dot_bracket(pdb_name:str, sequence:list, begin_end:list, output_dir:str, rev:bool=False, gnra_motif:str='GAAA'):
    """
    This function modifies the sequence by adding GNRA motif (e.g. GAAA) to join the two chains.
    Moreover, additional GC pair(s) are added to the sequence at the beginning and at the end.
    Together with the sequence, the 2D structure is also written to a file in dot-bracket format.
    The canonical pairings are given in the chainbraker positions (begin_end param).
    Args:
        pdb_name (str): pdb name
        sequence (list): list of sequences, e.g. ['CACGGCG', 'CGG']
        begin_end (list): list of chainbraker positions, e.g. [1, 0, 0, 0, 0, 0, 1, 1, 0, 1]
        output_dir (str): output directory
        rev (bool): if True, reverse the sequence
    """
    if rev:
        pdb_name = f'{pdb_name}_rev'
        begin_end = begin_end[::-1]
    
    # add GNRA motif to join the two chains
    seq = sequence[0] + gnra_motif + sequence[1]
    # add dots to the sequence
    str2d = ['.'] * len(seq)
    # add canonical pairings
    beg_end = np.where(np.array(begin_end) == 1)[0]  # get indeces where chainbraker is 1, e.g. [0, 6, 7, 9]
    beg_end_pairs = beg_end.reshape(-1, 2) # reshape to pairs, e.g. [[0, 6], [7, 9]]
    for i, pair in enumerate(beg_end_pairs):
        a, b = pair
        if i == 0: # strand 1
            str2d[a] = '('
            str2d[b] = '('
        elif i == 1: # strand 2
            str2d[a+len(gnra_motif)] = ')'
            str2d[b+len(gnra_motif)] = ')'
    str2d = ''.join(str2d)
    # add GC pair(s) to the beginning and the end of the sequence
    seq = 'GC' + seq + 'GC'
    str2d = '((' + str2d + '))'
    dot_file = f'>{pdb_name}\n{seq}\n{str2d}'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(f'{output_dir}/{pdb_name}.dot', 'w') as f:
        f.write(dot_file)

def write_commands(commands, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(f"{output_dir}/commands.txt", 'w') as f:
        for c in commands:
            f.write(f"{c}\n")

def get_centroids(path='data/il_3.78.json'):
    """
    Read json file and return first occurence as a centroid for each class
    """
    centroids = []
    with open(path) as f:
        json_file = json.load(f)
    for entry in json_file:
        for name in entry['alignment'].keys():
            centroids.append(name)
            break
    return centroids

def main():
    args = parse_args()
    if args.input is not None and args.directory is not None:
        print('Please specify either input file or directory')
        return
    
    centroids = get_centroids()
    if args.rna_tools is not None:
        print("Preparing commands for rna-tools")
    else:
        print("Preparing commands for RNAqua")
    
    if args.input is not None:
        cmds = process_file(args.input, args.output_dir, centroids=centroids, rna_tools = args.rna_tools)
    elif args.directory is not None:
        cmds = process_directory(args.directory, args.output_dir, centroids=centroids, rna_tools = args.rna_tools)
    else:
        print('Please specify either input file or directory')
        return
    assert args.output_dir is not None
    assert len(cmds) > 0
    print(f'Number of commands: {len(cmds)}')
    write_commands(cmds, args.output_dir)

if __name__ == '__main__':
    main()