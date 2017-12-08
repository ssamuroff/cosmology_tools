import numpy as np
import scipy, glob, argparse, yaml
import fitsio as fi
import pylab as plt
from tools.cosmosis import cosmosis_tools as ct


def main(args):
    ct.postprocess(args.chain, "output")
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--chain',"-c", type=str, action='store')

    args = parser.parse_args()

    main(args)
