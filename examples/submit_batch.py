import shutil
import os
from subprocess import check_output
from glob import glob
import time
import sys
sys.path.append('/n/home06/jkemp/SharedDMT')
from parse_io import read_args, write_args



def main(argv):
    """Submit a batch job for DMT.
- command: the program name
- args: the input parameter filename
- argname_to_vary: the argument name in the input parameter list to
                   vary between processes.
- npoint, begin, end: Calculate with  arg in a CLOSED interval [begin, end]
                      linearly spaced with npoint points.
- nprocess: The number of jobs (cores = jobs * cores per job) to ask the cluster for.
- [dry]: set to make see the batch_script but not submit.

If npoint is not a mutiple of nprocess, ADDS AN EXTRA JOB to deal
with the remainder (which will by definition have less work than the rest).
"""
    if len(argv) != 8 and len(argv) != 9:
        print("usage submit_batch.py command args argname_to_vary npoint begin end nprocess [dry]")
        sys.exit(2)

    command = argv[1]
    argfilename = argv[2]
    vary = argv[3]
    npoint = int(argv[4])
    begin = float(argv[5])
    end = float(argv[6])
    in_nproc = int(argv[7])
    nproc = in_nproc


    dry = True if len(argv) > 8 else False

    args = read_args(argfilename)
    outputname = args["OutputName"]

    if vary not in args:
        print("argname_to_vary not in argument file!")
        sys.exit(2)

    if nproc <= 1:
        print("Please do not use submit_batch for one processor (or less?)")
        sys.exit(2)

    if nproc > npoint:
        print("Warning: asked for more processes than points! "
              "Setting nprocess = npoint")
        nproc = npoint

    step = (end-begin)/(npoint-1)
    nperinterval = npoint//nproc
    endnperinterval = npoint-nproc*npoint
    if endnperinterval > 0:
        nproc = nproc +1
    pidstring = str(os.getpid())
    prefix = 'tmp/' + pidstring + '_' + argfilename + '_'
    for j in range(nproc):
        varyval = begin+j*nperinterval*step
        for k in range(nperinterval if j < in_nproc else endnperinterval):
            varyval = varyval + k*step
            args[vary] = str(varyval)
            args["OutputName"] = outputname + "_" + vary + "_" + str(varyval)
            write_args(args, prefix + str(j) + "_" + str(k))


    batch_options = """#!/bin/bash
#SBATCH -J dmt # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 1000 # Memory request (4Gb)
#SBATCH -t 0-0:10 # Maximum execution time (D-HH:MM)
#SBATCH -o tmp/dmt_%A_%a.out # Standard output
#SBATCH -e tmp/dmt_%A_%a.err # Standard error
module load gcc
module load intel-mkl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/n/sw/intel-cluster-studio-2019/compilers_and_libraries/linux/lib/intel64
"""
    batch_options += ("\nN=$(( ${SLURM_ARRAY_TASK_ID} == "+str(in_nproc+1)+" ? "+str(endnperinterval)+" : "+str(nperinterval)+" ))  " +
                      "\nfor ((i=0;i<N;i++)); do" +
                      "\n    " + command + " " + prefix + "${SLURM_ARRAY_TASK_ID}_${i}"
                      "\ndone"
    )
    print(batch_options)
    with open("batch_script",'w') as f:
        f.write(batch_options)
    if not dry:
        print(check_output(['sbatch','--array=0-'+str(nproc-1), "batch_script"], universal_newlines=True))
    print("Done.")


if __name__ == "__main__":
    main(sys.argv)
