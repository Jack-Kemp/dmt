import sys
sys.path.append('../')
import checkpoint_glob_all



def main(argv):
    """A redirect script to the main checkpoint_glob_all in SharedDMT
    folder.
    """
    checkpoint_glob_all.main(["",argv[1]])

if __name__ == "__main__":
    main(sys.argv)
