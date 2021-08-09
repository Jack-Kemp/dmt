import sys
import numpy as np
import glob
import checkpoint_glob


def main(argv):
    """Runs checkpoint_glob on all different sets of possible template
    file names of the form "template_name_t_*.dat. in a directory.
    """
    cdir = argv[1]
    if len(argv) != 2:
        print("usage checkpoint_glob_all.py dir_name.")
        sys.exit(2)
    datafs = glob.glob(cdir+"*_t_*.dat")
    datafs = [f.split("_t_")[0] for f in datafs]
    datafs = set(datafs)
    for f in datafs:
        print("Globbing: " + f)
        checkpoint_glob.main(["",f+"_t_*.dat"])

if __name__ == "__main__":
    main(sys.argv)
