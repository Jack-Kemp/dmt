import sys
from parse_io import ReadData, writeData
import numpy as np
import glob

def main(argv):
    if len(argv) != 2:
        print("usage slice.py name, will use first *")
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
    dfret = ReadData(datafs[-1])
    Tstart = dfret["t"][0]
    while Tstart > 0:
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
