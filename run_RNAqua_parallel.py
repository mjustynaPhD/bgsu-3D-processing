# Read RNAqua data and run parallel jobs
#
import os
import multiprocessing as mp
from tqdm import tqdm

CPUS = 4

def run_command(command):
    # command = f'cd ~/rnacomposer && {command}'
    # print(command)
    os.system(command)

def main():
    # read RNAqua data
    with open('RNAtools/commands.txt') as f:
        lines = f.readlines()
    
    # run commands (each line) in parallel
    # pool = mp.Pool(CPUS)
    # pool.map(run_command, lines)
    for l in tqdm(lines):
        run_command(l.strip())

if __name__ == '__main__':
    main()