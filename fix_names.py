import os
import sys

def main():
    args = sys.argv
    print(args)
    if len(args) < 2:
        print("Usage: fix_names.py <path>")
        sys.exit(1)
    path = args[1]

    # list all files in directory. If name contains '-' replace it with '_'
    files = os.listdir(path)
    for f in files:
        if '-' in f:
            new_name = f.replace('-', '_')
            os.rename(os.path.join(path, f), os.path.join(path, new_name))
            print(f'{f} -> {new_name}')
    print("Done")

if __name__ == '__main__':
    main()