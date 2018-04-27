#!/usr/bin/env python

# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob, argparse
import subprocess as sub

sys.path.append('/home/rpernak/python/GIT_python_modules')
import utils

# GLOBAL VARIABLES
DEFBUILDDIR = '/nas/ASR/LINEFILE_BUILD_TES_V_2.0'
DEFLNFL = \
  '/nas/project/rc/rc1/lnfl_local_version/lnfl_v3.1_linux_pgi_sgl'

class genTAPE3():
  """
  Assuming that line files exist for each HITRAN molecule, generate 
  an associated binary TAPE3
  """

  def __init__(self, lineFileDir=DEFBUILDDIR, \
    lnflPath=DEFLNFL, wnBounds=[0.0, 200000.0], \
    tape5Dir='%s/LNFL_input_files' % DEFBUILDDIR, \
    tape3Dir='%s/lnfl/TAPE3_files' % DEFBUILDDIR, \
    wvIso=False):
    """
    Constructor for genTAPE3 class -- where are ASCII line files, 
    what are the executable paths, where should output files go, etc.

    Inputs

    Outputs

    Keywords
      lineFileDir -- string, full path to line file build directory
        (where construction of a new line file is done 
         before releasing it); it is assumed that line_file/, 
         line_files_By_Molecule/, and lnfl/ subdirectories exist 
         underneath it
      lnflPath -- string, full path to LNFL executable to use
      wnBounds -- float array, 2-element array of starting and ending 
        wavenumbers for the TAPE3 files (25 cm-1 will be added to 
        each bound to include "relevant" line contributions)
      tape5Dir -- string, path to which TAPE5 files that will be 
        generated are written
      tape3Dir -- string, path to which TAPE3 files that will be 
        generated are written
      wvIso -- boolean, only process the water vapor isotopologues 
        (each of which has a separate line file)
    """

    lfMolDir = '%s/line_files_By_Molecule' % lineFileDir
    pathCheck = [lineFileDir, lfMolDir, lnflPath]
    for path in pathCheck: utils.file_check(path)

    # set some attributes necessary for the rest of the processing
    self.lineFileDir = str(lineFileDir)
    self.pathLNFL = str(lnflPath)
    self.mols = sorted(glob.glob('%s/*' % lfMolDir))
    self.nMols = len(self.mols)
    self.dirT5 = str(tape5Dir)
    self.dirT3 = str(tape3Dir)
    self.wnLims = list(wnBounds)
    self.isoH2O = bool(wvIso)

    # do the output dirs exist?
    for path in [tape5Dir, tape3Dir]: self.makeDirs(path)

  # end constructor

  def makeDirs(self, inPath):
    """
    Make a directory if it doesn't already exist
    """

    if not os.path.exists(inPath): os.mkdir(inPath)
  # end makeDirs()

  def cleanUp(self):
    """
    Remove any TAPE files that were generated
    """

    tapeList = sorted(glob.glob('TAPE?'))
    tapeList = ['TAPE%d' % num for num in [1, 2, 5, 6, 7, 10]]
    for tape in tapeList:
      if os.path.isfile(tape): os.remove(tape)
    # end TAPE loop

  # end cleanUp

  def makeLinks(self, source, target):
    """
    Make symbolic links

    Inputs
      source -- string, path to file to be linked
      target -- string, name associated with link
    """

    if os.path.exists(target): os.unlink(target)
    os.symlink(source, target)
  # end makeLinks()

  def makeTAPE5(self):
    """
    Make the TAPE5 LNFL specifications file for each molecule
    """

    wn1, wn2 = self.wnLims

    # loop through each HITRAN molecule and create an associated TAPE5
    allT5 = []
    for iMol, mol in enumerate(self.mols):
      base = os.path.basename(mol)
      print(base)
      tape5 = 'TAPE5_%s' % base

      # LNFL TAPE5 records 
      # (see lnfl_instructions document in LNFL release)
      rec1 = '$ %s' % base
      rec2 = '%10.3f%10.3f' % (wn1-25, wn2+25)

      # start off with all molecules off, then turn iMol on, then 
      # generate a single string instead of a list of characters
      # and append 
      rec3 = ['0'] * self.nMols
      rec3[iMol] = '1'
      rec3 = ''.join(rec3) + '     NBLK1 NOCPL LNOUT '
      end = '%%%%%'

      outDat = [rec1, rec2]

      # line coupling molecules
      if base in ['02_CO2', '06_CH4', '07_O2']:
        rec3 = rec3.replace('NOCPL', 'MRG2')
        rec4 = [' '] * self.nMols
        rec4[iMol] = '1'
        rec4 = ''.join(rec4)
        outDat += [rec3, rec4]
      else:
        outDat.append(rec3)
      # endif coupling

      outDat.append(end)

      # now write TAPE5
      outFP = open(tape5, 'w')
      for line in outDat: outFP.write('%s\n' % line)
      outFP.close()

      # copy TAPE5 to subdirectory for molecule in buildDir
      target = '%s/%s' % (self.dirT5, tape5)
      if os.path.exists(target):
        print('WARNING: overwriting %s' % target)
      # endif target check
      os.rename(tape5, target)

      allT5.append(target)
    # end molecule loop

    self.allT5 = list(allT5)
    return self

  # end makeTAPE5()

  def runLNFL(self):
    """
    Run LNFL for each molecule, thus generating a binary TAPE3
    """

    if 'allT5' not in dir(self):
      self.allT5 = sorted(glob.glob('%s/TAPE5_*' % self.dirT5))

    # set up the input directory
    self.makeLinks(self.pathLNFL, 'lnfl')
    tapeStrList = ['TAPE1', 'TAPE5']
    self.cleanUp()

    # loop through each HITRAN molecule and create an associated TAPE5
    for iMol, mol in enumerate(self.mols):
      base = os.path.basename(mol)
      print(base)
      tape5 = self.allT5[iMol]

      if self.isoH2O:
        # there are multiple line files to consider for H2O
        isoStr = ['01_h2o_161_only', '01_h2o_162_excl', \
          '01_h2o_162_only', '01_h2o_171_only', '01_h2o_172_only', \
          '01_h2o_181_only', '01_h2o_182_only']
        tape1List = ['%s/%s' % (mol, iso) for iso in isoStr]
      else:
        tape1List = ['%s/%s' % (mol, base)]
      # endif WV

      # loop really only exists for H2O
      for tape1 in tape1List:
        tapeList = [tape1, tape5]

        # grab the line coupling file if necessary
        if base in ['02_CO2', '06_CH4', '07_O2']:
          tape2 = '%s/lncpl_lines' % mol
          tapeList.append(tape2)
          tapeStrList.append('TAPE2')
        # endif line coupling

        # stage the files necessary for an LNFL run
        for source, target in zip(tapeList, tapeStrList):
          self.makeLinks(source, target)

        # call LNFL and save TAPE3 to unique name
        sub.call(['lnfl'])
        if self.isoH2O:
          tape3 = '%s/TAPE3_%s' % (mol, os.path.basename(tape1))
        else:
          tape3 = '%s/TAPE3_%s' % (mol, base)
        # endif wv
        if os.path.exists(tape3):
          print('WARNING: overwriting %s' % tape3)
        os.rename('TAPE3', tape3)

        # clean up
        self.cleanUp()
      # end TAPE1 loop

      #self.cleanUp()
      # if we're doing WV isotopologues, *only* do them
      if self.isoH2O: return
    # end molecule loop

    return
  # end runLNFL()
# end genTAPE3

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Generate a TAPE5 for each HITRAN molecule.')
  parser.add_argument('-dir', '--build_dir', type=str, \
    default=DEFBUILDDIR, \
    help='Line file build directory (which we assume contains a ' + \
    '"line_files_by_Molecule/" and "LNFL_input_files" subdirs.)')
  parser.add_argument('-lnfl', '--lnfl_path', type=str, \
    default=DEFLNFL, \
    help='Full path to LNFL executable to be used.')
  parser.add_argument('-t5', '--tape5_dir', type=str, \
    default='TAPE5_files', help='Directory for TAPE5 files.')
  parser.add_argument('-t3', '--tape3_dir', type=str, \
    default='TAPE3_files', \
    help='Relative path to directory for resulting TAPE3 files.')
  parser.add_argument('-wv', '--water_vapor', action='store_true', \
    help='Process the separate water vapor files (i.e., one for ' + \
    'each isotopologue that the JPL team wants).')
  parser.add_argument('-wn', '--wavenum_lims', nargs=2, \
    type=float, default=[500, 3500], \
    help='Wavenumber limits (cm-1) for TAPE3 (after LNFL run)')
  parser.add_argument('-5', '--only_tape5', action='store_true', \
    help='Only write the TAPE5s (no LNFL run).')
  parser.add_argument('-3', '--only_tape3', action='store_true', \
    help='Run LNFL without first making the TAPE5s.')
  args = parser.parse_args()

  tape3 = genTAPE3(lineFileDir=args.build_dir, \
    lnflPath=args.lnfl_path, wnBounds=args.wavenum_lims, \
    tape5Dir=args.tape5_dir, tape3Dir=args.tape3_dir, \
    wvIso=args.water_vapor)

  if args.only_tape5:
    tape3.makeTAPE5()
  elif args.only_tape3:
    tape3.runLNFL()
  else:
    tape3.makeTAPE5()
    tape3.runLNFL()
  # endif
# endif main()

