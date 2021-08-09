from collections import OrderedDict
from numpy import genfromtxt, savetxt
from numpy.lib.recfunctions import structured_to_unstructured


def _read_args_line(line, args):
    """Helper function for read_args. Reads a single line into string = string."""
    try:
        k, v = line.split('=')
        args[k.strip()] = v.strip()
    except:
        pass

def _parse_args_line(line, args):
    """Helper function for read_args. Parses a single line, that is,
    does its best to convert value to python variable rather than string.
    """
    try:
        k, v = line.split('=')
        v = v.strip()
        k = k.strip()
        try:
            args[k] = int(v)
        except ValueError:
            try:
                args[k] = float(v)
            except ValueError:
                if v == "true":
                    args[k] = True
                elif v == "false":
                    args[k] = False
                else:
                    args[k] = v
    except:
        pass

def read_args(name):
    """Read input parameters into OrderedDict."""
    with open(name, 'r') as f:
        args = OrderedDict()
        for line in f:
            _read_args_line(line, args)
    return args

def parse_args(name):
    """Parse input parameters into OrderedDict."""
    with open(name, 'r') as f:
        args = OrderedDict()
        for line in f:
            _parse_args_line(line, args)
    return args

def write_args(args, output):
    """Write input parameters from OrderedDict to output file output."""
    with open(output, 'w') as f:
        f.write("input\n{\n")
        for k, v in args.items():
            f.write(k + " = " +  str(v) + "\n")
        f.write("\n}")

def append_args(args, output, prefix='', input_header=True):
    """Write input parameters from OrderedDict to output file output."""
    with open(output, 'a') as f:
        if input_header:
            f.write(prefix+"input\n"+prefix+"{\n")
        for k, v in args.items():
            f.write(prefix+k + " = " +  str(v) + "\n")
        if input_header:
            f.write("\n"+prefix+"}")


class ReadData:
    """Read in data with ReadData(name). Gets data, input parameters and
    run information. Can be accessed like a dict based on variable name.
    For 2D data:
    -ReadData["name"] returns a matrix of all values.
    -ReadData["name", 3] returns column 3 (equivalent to ReadData["name_3"])
    -ReadData["name", 2:5] or ReadData["name", 2,3,4] returns matrix of columns 2,3,4.
    """
    def __init__(self, name, nSitesName='N'):
        self.data = genfromtxt(name, names=True)
        self.args = OrderedDict()
        self.runinfo = OrderedDict()
        with open(name, 'r') as f:
            for line in f:
                if line[0:2] == "#:":
                    _parse_args_line(line[2:], self.runinfo)
                elif line[0:2] == "#~":
                    _parse_args_line(line[2:], self.args)
        self.N = self.args[nSitesName]

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            try:
                return self.data[key]
            except ValueError:
                try:
                    return structured_to_unstructured(self.data[[key+"_"+str(i) for i in range(self.N)]])
                except KeyError:
                    try:
                        return self.args[key]
                    except KeyError:
                        return self.runinfo[key]
        else:
            if isinstance(key[-1], int):
                if len(key) == 2:
                    return self.data[key[0] + "_" + str(key[-1])]
                return structured_to_unstructured(self.data[[key[0]+ "_" + str(i) for i in key[1:]]])
            if isinstance(key[-1], slice):
                return structured_to_unstructured(self.data[[key[0]+ "_" + str(i) for i in range(*key[-1].indices(self.N))]])
            return self.data[key]


def writeData(name, rdata):
    cnames = rdata.data.dtype.names
    savetxt(name, rdata.data, header = ' '.join(cnames))
    append_args(rdata.runinfo, name, prefix='#:', input_header=False)
    append_args(rdata.args, name, prefix='#~')
