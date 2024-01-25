# Read RNAqua data and run parallel jobs
#
import os
import sys
import multiprocessing as mp

CPUS = 24

def run_command(command):
    command = f'cd ~/rnacomposer && {command}'
    os.system(command)

def main():
    # read RNAqua data
    with open('RNAqua') as f:
        lines = f.readlines()
    
    # run commands (each line) in parallel
    # pool = mp.Pool(CPUS)
    # pool.map(run_command, lines)
    for l in lines:
        run_command(l.strip())

if __name__ == '__main__':
    main()