import numpy as np
import scipy.integrate as integrate
import argparse
from ambcm_core import *


parser = argparse.ArgumentParser("ambcm")
parser.add_argument("z", help="Cos of angle", type=float)
parser.add_argument("lambd", help="Particle albedo", type=float)
parser.add_argument("type", help="simple, small, mid or auto, depending on your z, set at simple by default", type=str)
args = parser.parse_args()


if   (args.type=="simple"):
    print(ambrcm(args.z, args.lambd))
elif (args.type=="small"):
    print(ambrcm_small(args.z, args.lambd))
elif (args.type=="mid"):
    print(ambrcm_mid(args.z, args.lambd))
elif (args.type=="auto"):
    print(ambrcm_auto(args.z, args.lambd))
else:
    print(ambrcm(args.z, args.lambd)) 

