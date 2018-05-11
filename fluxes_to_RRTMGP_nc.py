#!/usr/bin/env python

from __future__ import print_function

import os, sys, glob, argparse
import numpy as np
import netCDF4 as nc

# this module was eventually moved to externals/common
#sys.path.append('externals/common')
import utils
import RC_utils as RC

"""
To do list:
  - remove need for an netCDF template
"""

class configSetup():
  def __init__(self, inFile):
    """
    Parse the input .ini file (inFile) and return as a dictionary for 
    use in the rest of this module

    Inputs
      inFile -- string, full path to .ini file that specifies ...

    Keywords
      doSW -- boolean, process shortwave instead of longwave
    """

    utils.file_check(inFile)

    # standard library, but name depends on Python version
    if sys.version_info.major < 3:
      import ConfigParser
    else:
      import configparser as ConfigParser
    # endif Python version

    cParse = ConfigParser.ConfigParser()
    cParse.read(inFile)
    cpSections = cParse.sections()

    # loop over each field (of all sections) and keep the field and 
    # associated value in returned object (self)
    for iCPS, cps in enumerate(cpSections):
      cItems = cParse.items(cps)
      for item in cItems: setattr(self, item[0], item[1])
    # end sections loop

  # end constructor
# end configSetup()

class swRRTMGP():
  def __init__(self, config, reverseVertical=False):
    """
    Class that conforms all of SW spectra generated by the SW Flux 
    Calculation software to the RRTMGP netCDF convention

    Input
      config -- configSetup object

    Keywords
    """

    paths = [config.top, config.nctemplate]
    for path in paths: utils.file_check(path)

    # gather all output netCDF files underneath the working 
    # directories generated in SW_create_inputs.py/SW_calc_fluxes.py
    search = '%s/%s/*.nc' % (config.top, config.subdirstr)
    profNC = sorted(glob.glob(search))
    self.profiles = profNC
    self.tempNC = config.nctemplate
    self.id = config.id
    self.specs = config.specs

    # for the new netCDF; we'll assume all profNC have same dimensions
    self.nProf = len(profNC)
    ncObj = nc.Dataset(profNC[0], 'r')
    self.nLev = ncObj.dimensions['levels'].size
    self.nLay = self.nLev - 1
    ncObj.close()

    # bands will remain unchanged
    ncObj = nc.Dataset(self.tempNC, 'r')
    self.bands = np.array(ncObj.variables['band_lims_wvn'])
    self.nBands = self.bands.shape[0]
    ncObj.close()
    
    self.base = os.path.basename(self.tempNC)

    subStr = 'sw'
    newSub = '%s-%s' % (subStr, self.id)
    if subStr in self.base:
      self.outFile = self.base.replace(subStr, newSub)
    else:
      self.outFile = '%s_%s' % (newSub, self.base)
    # endif base

    # fluxes are produced by LBLRTM, so we need to remove any RRTMGP
    # strings in the output filename
    self.outFile = self.outFile.replace('rrtmgp', 'lblrtm')

    # these three attribute lists have elements that correspond to 
    # each other ('down_direct' : 'band_flux_dir_dn' : 'flux_dir_dn')
    # but there's no reason to have RRTMGP wavenumbers because of 
    # self.bands
    self.fluxVars = ['down_direct', 'down_diffuse', \
      'down_total', 'up_total', 'net_flux', 'wavenumber']
    self.rrtmgpVars = ['band_flux_dir_dn', 'band_flux_dif_dn', \
      'band_flux_dn', 'band_flux_up', 'band_flux_net', '']
    self.rrtmgpVarsBB = ['flux_dir_dn', 'flux_dif_dn', \
      'flux_dn', 'flux_up', 'flux_net', '']

    self.sw = True
    self.reverse = reverseVertical
  # end constructor()

  def initializeNC(self):
    """
    Using the netCDF template, start a netCDF that will contain the 
    fields from the RRTMGP netCDF but will be stuffed with fluxes as 
    calculated with the SW Flux Calculation library

    the template is expected to follow the 
    rrtmgp-??-inputs-outputs-*.nc convention, where ?? is either "lw"
    or "sw"
    """

    print('Building %s' % self.outFile)

    inObj = nc.Dataset(self.tempNC, 'r')
    outObj = nc.Dataset(self.outFile, 'w')

    # for RRTMGP, we only have 5 dimensions, but lev, lay, and col
    # will probably be different in the SW flux calc output
    dims = inObj.dimensions
    for dim in dims:
      dimName = dims[dim].name
      if dimName == 'lev':
        dimSize = self.nLev
      elif dimName == 'lay':
        dimSize = self.nLay
      elif dimName == 'col':
        dimSize = self.nProf
      else:
        dimSize = dims[dim].size
      # endif dim.name

      outObj.createDimension(dim, dimSize)
    # end dim loop

    # unchanged, necessary variables from template
    inVar = inObj.variables['band_lims_wvn']
    outVar = outObj.createVariable(\
      inVar.name, inVar.dtype, inVar.dimensions)
    outVar.units = 'cm-1'
    outVar[:] = np.array(inVar)

    # there is no band dependence on P or TSI, so store them now
    inVarLev = inObj.variables['p_lev']
    inVarLay = inObj.variables['p_lay']
    inVarTSI = inObj.variables['total_solar_irradiance']
    if self.specs != '':
      specsObj = nc.Dataset(self.specs)
      pLev = np.array(specsObj.variables['pres_level'])
      pLay = np.array(specsObj.variables['pres_layer'])
      tsi = np.array(specsObj.variables['total_solar_irradiance'])
      pUnits = specsObj.variables['pres_level'].units
      tsiUnits = specsObj.variables['total_solar_irradiance'].units
      sza = np.array(specsObj.variables['solar_zenith_angle'])
      specsObj.close()

      if self.sw:
        # filter out nighttime profiles for SW
        # making some assumptions here that iUseProf.size and 
        # self.nProf are equal, and that there's a 1:1 correspondence
        iUseProf = np.where(sza < 90)[0]
        if iUseProf.size == 0: sys.exit('No daytime profiles found')
      else:
        iUseProf = np.arange(self.nProf)
      # endif SW

      outVarLev = outObj.createVariable(\
        inVarLev.name, inVarLev.dtype, inVarLev.dimensions)
      outVarLev.units = pUnits
      outVarLev[:] = pLev.T[::-1, iUseProf] if self.reverse else \
        pLev.T[:, iUseProf]

      outVarLay = outObj.createVariable(\
        inVarLay.name, inVarLay.dtype, inVarLay.dimensions)
      outVarLay.units = pUnits
      outVarLay[:] = pLay.T[::-1, iUseProf] if self.reverse else \
        pLay.T[:, iUseProf]

      outVarTSI = outObj.createVariable(\
        inVarTSI.name, inVarTSI.dtype, inVarTSI.dimensions)
      outVarTSI.units = tsiUnits
      outVarTSI[:] = tsi.T[iUseProf]
    else:
      # fill these in eventually
      pLev = np.zeros((self.nLev, self.nProf))
      pLay = np.zeros((self.nLay, self.nProf))
      pUnits = 'mbar'
    # endif specs

    # and now the variables that will change
    keys = self.fluxVars
    for iKey, key in enumerate(keys):
      if key == 'wavenumber': continue
      inVar = inObj.variables[self.rrtmgpVars[iKey]]
      inVarBB = inObj.variables[self.rrtmgpVarsBB[iKey]]

      outVar = outObj.createVariable(\
        inVar.name, inVar.dtype, inVar.dimensions)
      outVar.units = inVar.units
      outVarBB = outObj.createVariable(\
        inVarBB.name, inVarBB.dtype, inVarBB.dimensions)
      outVarBB.units = inVarBB.units
    # end keep loop

    # heating rates
    if self.sw:
      inVar = inObj.variables['band_heating_rate']
      outVar = outObj.createVariable(\
        inVar.name, inVar.dtype, inVar.dimensions)
      outVar.units = inVar.units

      inVarBB = inObj.variables['heating_rate']
      outVarBB = outObj.createVariable(\
        inVarBB.name, inVarBB.dtype, inVarBB.dimensions)
      outVarBB.units = inVarBB.units
    # endif sw

    inObj.close()
    outObj.close()

    return self
  # end initializeNC()

  def combineArr(self):
    """
    Combine flux arrays from SW Flux Calculation output files
    """

    # expected variable names
    keys = self.fluxVars

    # first make lists of flux and wavenumber arrays
    profDict = {}
    for key in keys: profDict[key] = []

    for prof in self.profiles:
      inObj = nc.Dataset(prof, 'r')

      for var in (inObj.variables):
        if var not in keys: continue
        profDict[var].append(np.array(inObj.variables[var]))
      # end var loop

      inObj.close()
    # end profile loop

    # wavenumber arrays should be identical
    ref = profDict['wavenumber'][0]
    for i in range(1, self.nProf):
      test = profDict['wavenumber'][i]
      errMsg = 'Wavenumbers for Experiments 1 and %d ' % (i+1) + \
        'are inconsistent, returning'
      if not np.all(ref == test): sys.exit(errMsg)
    # end profile loop

    # and since wavenumber arrays are identical (now), only need 1
    profDict['wavenumber'] = np.array(ref)

    # rearrange the flux arrays to match RRTMGP convenvtion of 
    # nLev x nProf x nBand and set appropriate attribute in object
    newDim = (1, 0, 2)
    for key in keys:
      newArr = np.array(profDict[key])
      if key == 'wavenumber':
        setattr(self, key, newArr)
      else:
        setattr(self, key, np.transpose(newArr, newDim))
      # endif key
    # end key loop
  # end combineArr()

  def computeBands(self, broadband=False, reverseVertical=False):
    """
    Compute fluxes for each RRTMGP-defined band

    Keywords
      broadband -- boolean, calculate broadband instead of by-band 
        fluxes
    """

    wnArr = self.wavenumber
    keys = self.fluxVars
    fluxDict = {}

    if broadband:
      for iKey, key in enumerate(keys):
        if key == 'wavenumber': continue
        fluxArr = getattr(self, key)
        fluxDict[key] = fluxArr[:, :, :].sum(axis=2)
        if self.reverse: fluxDict[key] = fluxDict[key][::-1, :]
      # end key loop

      rrtmgpVars = self.rrtmgpVarsBB
    else:
      for key in keys:
        # first initialize arrays in fluxDict
        if key == 'heat_rate':
          # only for LW, since HR is calculated by radsum
          fluxDict[key] = \
            np.zeros((self.nLay, self.nProf, self.nBands))
        else:
          fluxDict[key] = \
            np.zeros((self.nLev, self.nProf, self.nBands))
        # endif HR
      # end key loop

      # now stuff the arrays in fluxDict
      for iBand, band in enumerate(self.bands):
        # the SW flux calculation should go down to 100 cm-1
        # but the lowest wavenumber in RRTMGP is 820; we wanna include
        # all SW radiation, so we'll extend the first band as low as 
        # we can
        if iBand == 0:
          cond = wnArr < band[1]
        else:
          cond = (wnArr >= band[0]) & (wnArr < band[1])
        # endif band

        for key in keys:
          if key == 'wavenumber': continue
          fluxArr = getattr(self, key)
          fluxDict[key][:, :, iBand] = fluxArr[:, :, cond].sum(axis=2)
        # end key loop
      # end band loop

      for key in keys: 
        if self.reverse: fluxDict[key] = fluxDict[key][::-1, :, :]

      rrtmgpVars = self.rrtmgpVars

    # endif BB

    # now add arrays to netCDF
    outObj = nc.Dataset(self.outFile, 'r+')
    for iKey, key in enumerate(keys):
      if key == 'wavenumber': continue
      outVar = outObj.variables[rrtmgpVars[iKey]]
      outVar[:] = fluxDict[key]
    # end keep loop

    if self.sw:
      # heating rates calculations (not entirely sure if my method is 
      # correct...)
      if broadband:
        outVarBB = outObj.variables['heating_rate']
        diff = np.array(outObj.variables['flux_dn']) - \
          np.array(outObj.variables['flux_up'])
        hr = RC.fluxToHR(diff)
        outVarBB[:] = hr[::-1, :] if self.reverse else hr
        outObj.close()
      else:
        outVar = outObj.variables['band_heating_rate']
        diff = np.array(outObj.variables['band_flux_dn']) - \
          np.array(outObj.variables['band_flux_up'])
        hr = RC.fluxToHR(diff)
        outVar[:] = hr[::-1, :, :] if self.reverse else hr
      # endif BB
    # endif SW
  # end computeBands()
# end swRRTMGP()

class lwRRTMGP(swRRTMGP):
  def __init__(self, config):
    """
    Class that conforms all of LW spectra generated by RADSUM to the 
    RRTMGP netCDF convention

    Most of the work can be done with methods from the swRRTMGP class,
    but the inputs will be different, and thus the configuration file
    and associated config object will also differ. Only the 
    constructor and combineArr() methods should have to be redefined
    for LW

    Right now, this is not very flexible -- i'm working off what i 
    have from RFMIP

    Input
      config -- configSetup object

    Keywords

    """

    paths = [config.top, config.nctemplate]
    for path in paths: utils.file_check(path)

    # gather all RADSUM output files for both LW bands (10-2000 cm-1 
    # and 2000-3250 cm-1)
    searchB1 = '%s/10-2000/%s/LBL_Runs/%s_*/OUTPUT_RADSUM' % \
      (config.top, config.exp, config.subdirstr)
    searchB2 = '%s/2000-3250/%s/LBL_Runs/%s_*/OUTPUT_RADSUM' % \
      (config.top, config.exp, config.subdirstr)
    profB1 = sorted(glob.glob(searchB1))
    profB2 = sorted(glob.glob(searchB2))

    if len(profB1) != len(profB2): sys.exit('Returning: unequal ' +\
      'amounts of profiles for 10-2000 and 2000-3250 cm-1 bands')

    # i stupidly do not do padded string formatting with my column 
    # directories (RFMIP LW), so the first three columns are 1, 10, 
    # and 100, so i have to sort better
    colNum = [os.path.basename(os.path.dirname(prof)) \
      for prof in profB1]
    colNum = np.array([int(col.split('_')[1]) for col in colNum])
    iSort = np.argsort(colNum)

    self.profilesB1 = np.array(profB1)[iSort]
    self.profilesB2 = np.array(profB2)[iSort]
    self.tempNC = config.nctemplate
    self.id = config.id
    self.specs = config.specs

    # for the new netCDF; we'll assume all profiles have same dims
    # and that the same amount of profiles exist for both bands
    self.nProf = len(profB1)
    tempDict = RC.radsumRead(profB1[0])
    self.nLev = tempDict['down_flux'].shape[1]
    self.nLay = self.nLev - 1

    # bands will remain unchanged
    ncObj = nc.Dataset(self.tempNC, 'r')
    self.bands = np.array(ncObj.variables['band_lims_wvn'])
    self.nBands = self.bands.shape[0]
    ncObj.close()

    self.base = os.path.basename(self.tempNC)

    subStr = 'lw'
    newSub = '%s-%s' % (subStr, self.id)
    if subStr in self.base:
      self.outFile = self.base.replace(subStr, newSub)
    else:
      self.outFile = '%s_%s' % (newSub, self.base)
    # endif base

    # fluxes are produced by LBLRTM, so we need to remove any RRTMGP
    # strings in the output filename
    self.outFile = self.outFile.replace('rrtmgp', 'lblrtm')

    # these three attribute lists have elements that correspond to 
    # each other ('down_flux' : 'band_flux_dn' : 'flux_dn')
    # but there's no reason to have RRTMGP wavenumbers because of 
    # self.bands; also, for LW, heating rate is calculated in RADSUM
    self.fluxVars = ['up_flux', 'down_flux', \
      'net_flux', 'heat_rate', 'wavenumber']
    self.rrtmgpVars = ['band_flux_up', 'band_flux_dn', \
      'band_flux_net', 'band_heating_rate', '']
    self.rrtmgpVarsBB = ['flux_up', 'flux_dn', 'flux_net', \
      'heating_rate', '']

    self.sw = False
  # end constructor

  def combineArr(self, test=False):
    """
    Combine flux arrays from SW Flux Calculation output files
    """

    # expected variable names
    keys = self.fluxVars

    # first make lists of flux and wavenumber arrays
    profDict = {}
    for key in keys: profDict[key] = []

    # read in ASCII RADSUM output, combine 2 bands together into 
    # single arrays for each variable for a given profile
    # very time consuming at 1 cm-1 resolution
    # so for testing i'll just store and load into a .npz file
    npzFile = 'RFMIP_lw_fluxes.npz'
    if test:
      print('Loading %s' % npzFile)
      profDict = np.load(npzFile)['profDict'].item()
      #print(profDict['profDict'].item().keys())
    else:
      for pb1, pb2 in zip(self.profilesB1, self.profilesB2):
        print('Reading %s' % pb1)
        b1Dict = RC.radsumRead(pb1)
        print('Reading %s' % pb2)
        b2Dict = RC.radsumRead(pb2)

        for key in keys:
          if key == 'wavenumber': continue

          # stack a given array in the wavenumber dimension
          # i.e., join together 2 LW bands for given profile
          if key == 'heat_rate':
            # heat rates are assigned to layers, not levels, but 
            # RADSUM output is on levels, so there is a fill value
            # at the "TOA" for heating rate
            stack = np.vstack(\
              (b1Dict[key][:, 1:], b2Dict[key][:, 1:]))
          else:
            stack = np.vstack((b1Dict[key], b2Dict[key]))
          # endif key

          profDict[key].append(stack)
        # end key loop
      # end profile loop
      np.savez(npzFile.split('.')[0], profDict=profDict)
      print('Wrote %s' % npzFile)
      #print(profDict.keys())
    # endif test

    # wavenumber arrays should be identical, so we only assemble
    # a wavenumber array from one profile
    print('Assembling wavenumber array')
    b1WN = RC.radsumRead(self.profilesB1[0])['wavenumber1']
    b2WN = RC.radsumRead(self.profilesB2[0])['wavenumber1']
    profDict['wavenumber'] = np.hstack((b1WN, b2WN))

    # rearrange the flux arrays to match RRTMGP convenvtion of 
    # nLev x nProf x nBand and set appropriate attribute in object
    newDim = (2, 0, 1)
    for key in keys:
      newArr = np.array(profDict[key])
      if key == 'wavenumber':
        setattr(self, key, newArr)
      else:
        setattr(self, key, np.transpose(newArr, newDim))
      # endif key
    # end key loop
  # end combineArr()
# end lwRRTMGP

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Reorganize the output from the SW Flux ' + \
    'Calculation library -- which produces a subdirectory ' + \
    'for every processed profile and an associated netCDF ' + \
    'with flux spectra underneath the subdir -- into a single ' + \
    'netCDF that conforms to the RRTMGP convention so that ' + \
    'fluxes and heating rates for two models (e.g., LBLRTM vs. ' + \
    'RRTMGP) can be compared. NOTE: right now, the code is ' + \
    'optimized to handle RFMIP profiles and has not been tested ' + \
    'with other (e.g., Garand) specifications.')
  parser.add_argument('--ini_file', '-ini', type=str, \
    default='sw_make_RRTMGP_nc.ini', \
    help='Configuration file with specs for ncRRTMGP class.' )
  parser.add_argument('--longwave', '-lw', action='store_true', \
    help='Process longwave output, which likely used radsum ' + \
    'output instead of output from the SW Flux Calculation repo.')
  parser.add_argument('--lw_test', '-test', action='store_true', \
    help='LW testing (involves some .npz saving and loading')
  args = parser.parse_args()

  iniFile = args.ini_file
  ini = configSetup(iniFile)

  if args.longwave:
    ncRRTMGP = lwRRTMGP(ini, reverseVertical=True)
    ncRRTMGP.initializeNC()
    ncRRTMGP.combineArr(test=args.lw_test)
    ncRRTMGP.computeBands()
    ncRRTMGP.computeBands(broadband=True)
  else:
    ncRRTMGP = swRRTMGP(ini, reverseVertical=True)
    ncRRTMGP.initializeNC()
    ncRRTMGP.combineArr()
    ncRRTMGP.computeBands()
    ncRRTMGP.computeBands(broadband=True)
  # endif LW
# end main()

