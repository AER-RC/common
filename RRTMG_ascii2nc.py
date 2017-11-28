#!/usr/bin/env python

import os, sys, glob, argparse
import numpy as np
import netCDF4 as nc
from nco import Nco as NCO

# RC GitLab repo
sys.path.append('common_modules')
import utils
import RC_utils as rc

# global variables
MODES = ['garand', 'rfmip']

class ascii():
  def findProfiles(self):
    """
    Extract the RRTMG flux files for each spectral domain (LW and SW)
    """

    lwFiles = sorted(glob.glob('%s/%s*' % (self.lwDir, self.search) ))
    swFiles = sorted(glob.glob('%s/%s*' % (self.swDir, self.search) ))
    if len(lwFiles) != len(swFiles):
      errMsg = 'LW and SW output were generated for different ' + \
        'amounts of profiles, returning'
      sys.exit(errMsg)
    # endif len

    # by now, if the number of LW files is zero, then so is the 
    # number of SW files
    if len(lwFiles) == 0: sys.exit('No profiles found, returning')

    self.nProfiles = len(lwFiles)

    return {'lw_files': lwFiles, 'sw_files': swFiles}
  # end findProfiles()

  def readASCII(self, inFile, shortWave=False):
    """
    Read a single RRTMG ASCII flux file and return the model output 
    in a dictionary to be used in makeNC()

    Input
      inFile -- string, full path to RRTMG ASCII flux file for a 
        single profile

    Output
      outDict -- dictionary with the following fields:
        level_pressure: pressure at layer boundaries (nLevel array)
        up_flux: upwelling flux (W/m2) as a function of wavenumber 
          and level (nLevel x nWavenumber array)
        net_flux: net flux (W/m2) as a function of wavenumber and
          level (nLevel x nWavenumber array)
        heat_rate: heating rate (K/day) as a function of wavenumber 
          and level (nLevel x nWavenumber array)
        wavenumber: spectral points (cm-1) vector (2 x nWavenumber)
          ([starting wavenumber of band, ending wavenumber of band])

        Longwave and Shortwave:
          down_flux: total downwelling flux (W/m2) as a function of 
            wavenumber and level (nLevel x nWavenumber array)

        Shortwave only:
          difdown_flux: diffuse down flux (W/m2) as a function of 
            wavenumber and level (nLevel x nWavenumber array)
          dirdown_flux: direct down flux (W/m2) as a function of 
            wavenumber and level (nLevel x nWavenumber array)

        NOTE: each flux field also has an associated broadband key 
          that contains an nLevel-element array of fluxes integrated 
          over the entire spectral domain

    Keywords
      shortWave -- boolean, process SW flux files instead of LW
    """

    profDat = open(inFile).read().splitlines()

    # these lists will include all spectral points and the broadband
    # pressure should be the same regardless of band
    pLev, wn1, wn2, upTot, downTot, net, hr, downDir, downDif = \
      ([] for i in range(9))

    # these lists are for single bands
    upTotBand, downTotBand, netBand, hrBand, dirBand, difBand = \
      ([] for i in range(6))
    pLevBand = []

    for line in profDat:
      split = line.split()

      if len(split) == 0:
        # empty lines imply the end of one band and the start of 
        # another, so we need to save the lists of fluxes (and HR) 
        # from the previous band and reset them for the new one
        if len(upTotBand) != 0:
          upTot.append(upTotBand)
          downTot.append(downTotBand)
          net.append(netBand)
          hr.append(hrBand)
          downDir.append(dirBand)
          downDif.append(difBand)

          # pLev can just be overwritten every band -- 
          # they are constant
          pLev = list(pLevBand)
        # end len check

        upTotBand, downTotBand, netBand, hrBand, dirBand, difBand = \
          ([] for i in range(6))

        # pLev can just be overwritten every band -- they are constant
        pLevBand = []

      elif split[0] == 'mb':
        # can skip this header info
        continue
      elif split[0] == 'LEVEL':
        # can skip this header info
        continue
      elif split[0] == 'Wavenumbers:':
        # extract spectral range of band then move to next line
        wn1.append(float(split[1]))
        wn2.append(float(split[3]))
        continue
      elif split[0] == 'Modules':
        # this is the footer, and we do not need anything from it
        break
      else:
        split = [float(i) for i in split]
        pLevBand.append(split[1])
        if len(split) == 6:
          # LW output
          upTotBand.append(split[2])
          downTotBand.append(split[3])
          netBand.append(split[4])
          hrBand.append(split[5])
        elif len(split) == 8:
          # SW output
          upTotBand.append(split[2])
          difBand.append(split[3])
          dirBand.append(split[4])
          downTotBand.append(split[5])
          netBand.append(split[6])
          hrBand.append(split[7])
        # end LW/SW

      # endif wn construction

    # end inDat loop

    outDict = {}

    outDict['wavenumber'] = np.array([wn1, wn2])[:, 1:].T
    outDict['level_pressures'] = np.array(pLev)

    # transpose the output arrays to follow RRTMGP netCDF convention
    # and slice to separate broadband from band arrays
    outDict['up_flux'] = np.array(upTot)[1:, :].T
    outDict['up_flux_BB'] = np.array(upTot)[0].T
    outDict['net_flux'] = np.array(net)[1:, :].T
    outDict['net_flux_BB'] = np.array(net)[0].T
    outDict['heat_rate'] = np.array(hr)[1:, :].T
    outDict['heat_rate_BB'] = np.array(hr)[0].T
    outDict['down_flux'] = np.array(downTot)[1:, :].T
    outDict['down_flux_BB'] = np.array(downTot)[0].T

    if shortWave:
      outDict['difdown_flux'] = np.array(downDif)[1:, :].T
      outDict['difdown_flux_BB'] = np.array(downDif)[0].T
      outDict['dirdown_flux'] = np.array(downDir)[1:, :].T
      outDict['dirdown_flux_BB'] = np.array(downDir)[0].T
    # end SW

    # spectrally sort the band data (not broadband)
    iSort = np.argsort(outDict['wavenumber'][:, 0])
    for key in outDict.keys():
      if 'BB' in key: continue
      if key in ['wavenumber', 'level_pressures']: continue
      outDict[key] = outDict[key][iSort, :]
    # end key loop
    outDict['wavenumber'] = outDict['wavenumber'][iSort, :]

    return outDict
  # end readASCII()

  def __init__(self, inDirLW, inDirSW, searchStr='OUTPUT_RRTM'):
    """
    Extract the RRTMG flux files for each spectral domain (LW and SW)

    Read a single RRTMG ASCII flux file and return the model output 
    in a dictionary to be used in makeNC()

    Input
      inDirLW -- string, directory with LW RRTMG files
      inDirSW -- string, directory with SW RRTMG files

    Output
      outDict -- dictionary with the following keys:
        lw_files: RRTMG LW ASCII file for each profile
        sw_files: RRTMG SW ASCII file for each profile

    Keywords
      searchStr -- string used for finding RRTMG ASCII files
    """

    self.lwDir = inDirLW
    self.swDir = inDirSW
    self.search = searchStr
    self.txtFiles = self.findProfiles()

    # read the LW ASCII files and store each in comprehensive dict
    lwDict = {}
    for iProf, prof in enumerate(self.txtFiles['lw_files']):
      lwDict['profile%03d' % (iProf+1)] = self.readASCII(prof)

    # read the SW ASCII files and store each in comprehensive dict
    swDict = {}
    for iProf, prof in enumerate(self.txtFiles['sw_files']):
      swDict['profile%03d' % (iProf+1)] = \
        self.readASCII(prof, shortWave=True)
    # end SW Loop

    self.lwFluxes = dict(lwDict)
    self.swFluxes = dict(swDict)

    # we now assume that all profiles have the same number of levels
    # and that the number is the same for each both spectral domains
    self.nLevels = lwDict['profile001']['level_pressures'].shape[0]
    self.nLayers = self.nLevels - 1
    self.pLev = lwDict['profile001']['level_pressures']

    self.lwWN = lwDict['profile001']['wavenumber']
    self.swWN = swDict['profile001']['wavenumber']
  # end constructor
# end ascii()

def writeNC(inObj, mode='garand', suffix='inputs-outputs.nc'):
  """
  Write a netCDF with the data in an 'ascii' object. This is done for
  each spectral domain (lw and sw)

  Input
    inObj -- object from ascii class

  Output
    None returned. netCDF is written to 'rrtmg-sw-suffix' and 
      'rrtmg-lw-suffix'

  Keywords
    mode -- string that determines what netCDF format to use
    suffix -- output netCDF filename suffix appended to 'rrtmg-?w-'
  """
  
# end writeNC()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Convert ASCII RRTMG output to netCDF format.')
  parser.add_argument('lw_dir', type=str, \
    help='Directory with RRTGM LW results.')
  parser.add_argument('sw_dir', type=str, \
    help='Directory with RRTGM LW results.')
  parser.add_argument('mode', type=str, \
    help='String that directs the script on what netCDF format ' + \
      'to use.')
  parser.add_argument('--search', type=str, default='OUTPUT_RRTM', \
    help='Search string that will be used to find RRTMG output ' + \
    'ASCII files.')
  args = parser.parse_args()

  lwDir = args.lw_dir; utils.file_check(lwDir)
  swDir = args.sw_dir; utils.file_check(swDir)

  ncMode = args.mode.lower()
  if ncMode not in MODES: sys.exit('Set mode to any of %s' % MODES)

  asciiObj = ascii(lwDir, swDir, searchStr=args.search)

  # did sanity checks for all fluxes and HR in Garand 1, all bands
  # LW sanity check
  #print(asciiObj.lwFluxes['profile001']['net_flux_BB'])
  """
  print(asciiObj.lwFluxes['profile001']['net_flux'][0])
  print()
  print(asciiObj.lwFluxes['profile001']['net_flux'][-1])
  """

  # SW sanity check
  #print(asciiObj.swFluxes['profile001']['heat_rate_BB'])
  """
  print(asciiObj.swFluxes['profile001']['heat_rate'][0])
  print()
  print(asciiObj.swFluxes['profile001']['heat_rate'][-1])
  """

# end main()

