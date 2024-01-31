# Read RNAqua data and run parallel jobs
#
import os
import multiprocessing as mp
from tqdm import tqdm

CPUS = 4
chains_errors = [
    "4V88",
    "4WSM",
    "6CZR",
    "7QR3",
    "7QR4",
    "7SZU",
    "7VTI",
    "7VYX",
    "7WIE",
    "8AF0",
    "8DK7",
    "8F4O",
    "8GXB",
    "8HBA",
    "8S95"
]

def run_command(command):
    # command = f'cd ~/rnacomposer && {command}'
    # print(command)
    os.system(command)

def main():
    # read RNAqua data
    with open('RNAtools/commands.txt') as f:
        lines = f.readlines()
    
    lines = [l.split(' ')[11][1:5] for l in lines]
    lines = [l for l in lines if l in chains_errors]

    # run commands (each line) in parallel
    # pool = mp.Pool(CPUS)
    # pool.map(run_command, lines)
    for l in tqdm(lines):
        run_command(l.strip())

if __name__ == '__main__':
    main()