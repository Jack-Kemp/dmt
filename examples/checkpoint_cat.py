import sys
sys.path.append('../')
import checkpoint_glob



def main(argv):
    """A redirect script to the main checkpoint_glob script in SharedDMT
    folder.
    """
    checkpoint_glob.main(["",argv[1]])

if __name__ == "__main__":
    main(sys.argv)
