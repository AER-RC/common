#!/usr/bin/env python

import os, sys, glob, argparse
import subprocess as sub
import numpy as np
import netCDF4 as nc

# conda install -c conda-forge pynco
from nco import Nco as NCO

# RC GitLab repo
# git clone git@lex-gitlab.aer.com:RC/common_modules.git
sys.path.append('common_modules')
import utils

# global variables (in CAPS)
MODES = ['garand', 'rfmip']

# LW and SW fields from RRTMGP (Garand) template that will be modified 
NCFIELDSG = ['band_heating_rate', 'band_lims_wvn', 'p_lev', \
  'band_flux_dn', 'band_flux_net', 'band_flux_up', \
  'flux_dn', 'flux_net', 'flux_up', 'heating_rate']

# SW fields from RRTMGP template that will be modified 
NCSWFIELDSG = ['band_flux_dif_dn', 'band_flux_dir_dn', \
  'flux_dif_dn', 'flux_dir_dn']

# LW and SW fields from RFMIP template that will be modified 
NCFIELDSR = []

# SW fields from RFMIP template that will be modified 
NCSWFIELDSR = []

# Directory that contains the RRTMGP netCDFs that will be used as a 
# template for the Garand atmospheres
REFNCDIR = '/rd47/scratch/RRTMGP/RRTMGP_SVN/trunk/' + \
  'test/flux_compute/ref'
REFNCLW = '%s/rrtmgp-lw-inputs-outputs-clear.nc' % REFNCDIR
REFNCSW = '%s/rrtmgp-sw-inputs-outputs-clear.nc' % REFNCDIR

class rrtmg():
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

    # some RRTMG files are stored with profile numbers that are not 
    # 0-padded (e.g., OUTPUT_RRTM.GARAND_1 instead of 
    # OUTPUT_RRTM.GARAND_01), so we need to try and address this
    if self.profiles == 'garand':
      fileList = [lwFiles, swFiles]
      for iList, fList in enumerate(fileList):
        profNum = np.array(\
          [int(prof.split('_')[-1]) for prof in fList])
        iSort = np.argsort(profNum)
        if iList == 0:
          lwFiles = list(np.array(lwFiles)[iSort])
        else:
          swFiles = list(np.array(swFiles)[iSort])
      # end fList loop
      
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
      outDict[key] = outDict[key][:, iSort]
    # end key loop
    outDict['wavenumber'] = outDict['wavenumber'][iSort, :]

    return outDict
  # end readASCII()

  def combineProfiles(self):
    """
    Merge together the fluxes and heating rates from all profiles 
    into a single nLevel x nProfile x nBand array for each parameter
    """

    fluxDict = [self.fluxesLW, self.fluxesSW]

    # loop over spectral domains and combine profiles for each
    for iFD, fDict in enumerate(fluxDict):
      pLev = []
      upTot, downTot, net, hr, downDir, downDif = \
        ([] for i in range(6))
      upTotBB, downTotBB, netBB, hrBB, downDirBB, downDifBB = \
        ([] for i in range(6))

      # loop over profiles
      for iKey, key in enumerate(fDict.keys()):
        pLev.append(fDict[key]['level_pressures'])

        # by-band fluxes (LW and SW)
        upTot.append(fDict[key]['up_flux'])
        downTot.append(fDict[key]['down_flux'])
        net.append(fDict[key]['net_flux'])
        hr.append(fDict[key]['heat_rate'])

        # broadband (LW and SW)
        upTotBB.append(fDict[key]['up_flux_BB'])
        downTotBB.append(fDict[key]['down_flux_BB'])
        netBB.append(fDict[key]['net_flux_BB'])
        hrBB.append(fDict[key]['heat_rate_BB'])

        # SW by-band and broadband
        if iFD == 1:
          downDir.append(fDict[key]['dirdown_flux'])
          downDif.append(fDict[key]['difdown_flux'])
          downDirBB.append(fDict[key]['dirdown_flux_BB'])
          downDifBB.append(fDict[key]['difdown_flux_BB'])
        # endif SW
      # end fDict keys loop

      # for transforming arrays from (nProfiles x nLevels x nBands)
      # to (nLevels x nProfiles x nBands)
      tAxes = (1,0,2)

      # now convert to arrays and assign as attributes to object
      combined = {}
      fields = self.ncFieldsSW if iFD == 1 else self.ncFieldsLW

      combined[fields[0]] = np.transpose(np.array(hr), axes=tAxes)
      combined[fields[2]] = np.array(pLev).T
      combined[fields[3]] = np.transpose(np.array(downTot), \
        axes=tAxes)
      combined[fields[4]] = np.transpose(np.array(net), axes=tAxes)
      combined[fields[5]] = np.transpose(np.array(upTot), axes=tAxes)
      combined[fields[6]] = np.array(downTotBB).T
      combined[fields[7]] = np.array(netBB).T
      combined[fields[8]] = np.array(upTotBB).T
      combined[fields[9]] = np.array(hrBB).T

      if iFD == 1:
        combined[fields[1]] = np.array(self.wnSW)
        combined[fields[10]] = np.transpose(np.array(downDif), \
          axes=tAxes)
        combined[fields[11]] = np.transpose(np.array(downDir), \
          axes=tAxes)
        combined[fields[12]] = np.array(downDifBB).T
        combined[fields[13]] = np.array(downDirBB).T
        self.combinedSW = dict(combined)
      else:
        combined[fields[1]] = np.array(self.wnLW)
        self.combinedLW = dict(combined)
      # endif SW/LW
    # end loop over spectral domains

    return self
  # end combineProfiles()

  def writeNC(self):
    """
    Write a netCDF with the data in an rrtmg object. This is done for
    each spectral domain (lw and sw)

      None returned. netCDF is written to 'rrtmg-sw-suffix' and 
        'rrtmg-lw-suffix'

    """

    # first copy over the netCDF templates
    lwCmd = [self.ncCopy, self.ncTempLW, self.ncOutLW]
    swCmd = [self.ncCopy, self.ncTempSW, self.ncOutSW]
    sub.call(lwCmd)
    sub.call(swCmd)

    # now edit the copies with profile data
    #print(self.combinedSW.keys())

    # for now, we'll ignore heating rate because the RRTMG output is 
    # on levels and the RRTMGP output is on layers
    tempFields = self.ncFieldsLW
    fields = list(tempFields)
    fields.remove(tempFields[0])
    fields.remove(tempFields[9])

    lwObj = nc.Dataset(self.ncOutLW, 'r+')
    swObj = nc.Dataset(self.ncOutSW, 'r+')

    print(np.array(lwObj.variables['flux_up']))
    print()
    print(self.combinedLW['flux_up'])
    #for field in fields:
    #  print(field)
    #  lwObj.variables[field][:] == self.combinedLW[field]
    # end fields loop

    lwObj.close()
    swObj.close()
  # end writeNC()

  def __init__(self, inDirLW, inDirSW, searchStr='OUTPUT_RRTM', \
    profiles='garand', templateLW=REFNCLW, templateSW=REFNCSW, \
    suffix='inputs-outputs.nc', ncCopyPath='nccopy'):
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
      profiles -- string that dictates what netCDF format is used 
        (e.g., Garand, RFMIP, etc.)
      templateLW -- string, full path to LW netCDF file for the 
        specified profiles
      templateSW -- string, full path to SW netCDF file for the 
        specified profiles
      suffix -- string, that is appended to "rrtmg-lw" and "rrtmg-sw" 
        in the output netCDF files
      ncCopyPath -- string, full path to nccopy executable (or just 
        "nccopy" if it is in $PATH)
    """

    self.lwDir = inDirLW
    self.swDir = inDirSW
    self.search = searchStr
    self.profiles = profiles
    self.txtFiles = self.findProfiles()
    self.ncTempLW = templateLW
    self.ncTempSW = templateSW
    self.ncCopy = ncCopyPath

    # determine which netCDF fields to modify
    if profiles == 'garand':
      self.ncFieldsLW = list(NCFIELDSG)
      self.ncFieldsSW = NCFIELDSG + NCSWFIELDSG
    else:
      self.ncFieldsLW = list(NCFIELDSR)
      self.ncFieldsSW = NCFIELDSR + NCSWFIELDSR
    # endif profiles

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

    self.fluxesLW = dict(lwDict)
    self.fluxesSW = dict(swDict)

    # we now assume that all profiles have the same number of levels
    # and that the number is the same for each both spectral domains
    self.nLevels = lwDict['profile001']['level_pressures'].shape[0]
    self.nLayers = self.nLevels - 1
    self.pLev = lwDict['profile001']['level_pressures']

    self.wnLW = lwDict['profile001']['wavenumber']
    self.wnSW = swDict['profile001']['wavenumber']
    self.nBandsLW = lwDict['profile001']['wavenumber'].shape[0]
    self.nBandsSW = swDict['profile001']['wavenumber'].shape[0]

    # now merge the profiles
    self.combineProfiles()

    # output filename construction
    self.ncOutLW = 'rrtmg-lw-%s' % suffix
    self.ncOutSW = 'rrtmg-sw-%s' % suffix

  # end constructor
# end rrtmg()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Convert ASCII RRTMG output to netCDF format.  ' + \
    'A netCDF for both the LW and SW is written to working directory.')
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
  parser.add_argument('--suffix', type=str, \
    default='inputs-outputs.nc', \
    help='Output netCDF filename suffix appended to "rrtmg-?w-" ' + \
    '(so for the default "inputs-outputs.nc" and for the LW, the ' + \
    'output netCDF filename would be "rrtmg-lw-inputs-outputs.nc.")')
  parser.add_argument('--lw_template', type=str, default=REFNCLW, \
    help='Full path to netCDF that will be used as a template ' + \
    'on which the output LW netCDF will be based.')
  parser.add_argument('--sw_template', type=str, default=REFNCSW, \
    help='Full path to netCDF that will be used as a template ' + \
    'on which the output SW netCDF will be based.')
  parser.add_argument('-n', '--nccopy_path', type=str, \
    default='/nas/project/p1770/dependencies/bin/nccopy', \
    help='Full path to the nccopy executable in the C netCDF ' + \
    'library (must be version 4.3.0 or newer).')
  args = parser.parse_args()

  lwDir = args.lw_dir; utils.file_check(lwDir)
  swDir = args.sw_dir; utils.file_check(swDir)

  ncMode = args.mode.lower()
  if ncMode not in MODES: sys.exit('Set mode to any of %s' % MODES)

  ncTempLW = args.lw_template; utils.file_check(ncTempLW)
  ncTempSW = args.sw_template; utils.file_check(ncTempSW)
  ncCopy = args.nccopy_path; utils.file_check(ncCopy)

  rrtmgObj = rrtmg(lwDir, swDir, searchStr=args.search, \
    profiles=ncMode, templateLW=ncTempLW, templateSW=ncTempSW, \
    suffix='inputs-outputs.nc', ncCopyPath=ncCopy)
  rrtmgObj.writeNC()

  # did sanity checks for all fluxes and HR in Garand 1, all bands
  # LW sanity check
  #print(rrtmgObj.fluxesLW['profile001']['net_flux_BB'])
  """
  print(rrtmgObj.fluxesLW['profile001']['net_flux'][0])
  print()
  print(rrtmgObj.fluxesLW['profile001']['net_flux'][-1])
  """

  # SW sanity check
  #print(rrtmgObj.fluxesSW['profile001']['heat_rate_BB'])
  """
  print(rrtmgObj.fluxesSW['profile001']['heat_rate'][0])
  print()
  print(rrtmgObj.fluxesSW['profile001']['heat_rate'][-1])
  """

# end main()

