import scm.plams as plams
import struct_generator2, paths, os, utility, job_results3, time, sys, geometry
from utility import bohr2angstrom
b2a = bohr2angstrom(1)

join = os.path.join