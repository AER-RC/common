#!/usr/bin/env python

import argparse
import utils
import FortranFile
from lblTools import readTape12

parser = argparse.ArgumentParser(\
  description='Quick test of Scott Zaccheo readTape12() with ' + \
  'single precision TAPE11 and a double precision TAPE12.')
parser.add_argument('-t11', '--tape11', action='store_true')
parser.add_argument('-t12', '--tape12', action='store_true')
args = parser.parse_args()

t11 = '/rd47/scratch/RRTMGP/band_generation/SW/diffuse_flux_calc/fluxes/TAPE11_direct_down_ndm_032'
t12 = '/nas/ASR/RC_Release_Benchmark_Tests/Solar/upwelling/output_files/TAPE12_LBLRTM_v12.8_AER_v3.6'

if args.tape11:
  outWN, param = readTape12(t11, double=False)
  print utils.pmm(param)
elif args.tape12:
  outWN, param = readTape12(t12, double=True)
  print utils.pmm(outWN)
  print utils.pmm(param)
  import matplotlib.pyplot as plot
  plot.plot(outWN, param, 'r')
  plot.savefig('t12_debug_readBinary.png')
else:
  print 'Did nothing, set --tape11 or --tape12'

