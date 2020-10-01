import sys
import shutil
import os
from subprocess import check_output
from glob import glob
import time

def main(argv):
    if len(argv) != 3:
        print("usage submit.py command args")
        sys.exit(2)
    command = argv[1]
    args = argv[2]
    pidstring = str(os.getpid())
    label = pidstring + '_'
    tmpargs = 'tmp/' + label + args
    shutil.copyfile(args,  tmpargs)

    batch_options = """#!/bin/bash
#SBATCH -J dmt # A single job name for the array
#SBATCH -n 4 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 8000 # Memory request (4Gb)
#SBATCH -t 0-18:00 # Maximum execution time (D-HH:MM)
#SBATCH -o tmp/dmt_%A_%a.out # Standard output
#SBATCH -e tmp/dmt_%A_%a.err # Standard error
module load gcc
module load intel-mkl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/n/sw/intel-cluster-studio-2019/compilers_and_libraries/linux/lib/intel64
"""
    batch_options += (
                      "\n" + command + " " + tmpargs
    )
    print(batch_options)
    with open("batch_script",'w') as f:
        f.write(batch_options)
    print(check_output(['sbatch','--array=0-'+str(1-1), "batch_script"], universal_newlines=True))
    print("Done.")


if __name__ == "__main__":
    main(sys.argv)
