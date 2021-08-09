import sys
from parse_io import ReadData, writeData
import numpy as np
import glob

def main(argv):
    """Takes a template filename of the form "checkpoint_t_*.dat", which
    represents a set of checkpoint output data from DMT, and combines all
    the unique data into one data file with continous time series, in case
    of restarts. The output file is "template_file_name + _globed".
    """
    if len(argv) != 2:
        print("usage check_point_glob.py name, will use first * in file name.")
        sys.exit(2)
    name = argv[1]
    datafs = glob.glob(name)
    datsuffix = None
    if name[-4:] == ".dat":
        name = name[:-4]
        datsuffix = -4
    Tind = name.split('_').index('*')
    T = np.array([float(f.split('_')[Tind][:datsuffix]) for f in datafs])
    T, datafs = (list(x) for x in zip(*sorted(zip(T, datafs))))
    ind = len(datafs)-1
    dfret = ReadData(datafs[ind])
    Tstart = dfret["t"][0]
    while Tstart > 0 and ind > 0:
        print(Tstart)
        ind = np.searchsorted(T, Tstart)
        print(ind)
        print(T[ind])
        dftest = ReadData(datafs[ind])
        dftest.data = np.delete(dftest.data,-1)
        dfret.data = np.hstack((dftest.data, dfret.data))
        Tstart = dfret["t"][0]
    writeData(name + '_globbed', dfret)

if __name__ == "__main__":
    main(sys.argv)
