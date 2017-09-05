import os, sys
import numpy as np

def readTAPE28(inFile, nSkip=52):
  """
  Read in a single TAPE28 (LBLRTM Brightness Temperature output file)
  and return spectrum in dictionary

  Call
    outDict = readTAPE28(inFile)

  Input
    inFile -- string, path to TAPE28 file

  Output
    outDict -- dictionary with wavenumber and brightness_temperature
      key/value pairs

  Keywords
    nSkip -- int, the number of lines in the header
  """

  # very, very simple function -- maybe just easier to use np.loadtxt
  # directly
  waveNum, bt = np.loadtxt(inFile, unpack=True, skiprows=nSkip)
  outDict = \
    {'wavenumber': waveNum , 'brightness_temperature': bt, 'units': 'K'}

  return outDict
# end readTAPE28

def readTAPE27(inFile, nSkip=52):
  """
  Read in a single TAPE27 (LBLRTM radiance output file)
  and return spectrum in dictionary

  Call
    outDict = readTAPE27(inFile)

  Input
    inFile -- string, path to TAPE28 file

  Output
    outDict -- dictionary with wavenumber and brightness_temperature
      key/value pairs

  Keywords
    nSkip -- int, the number of lines in the header
  """

  # very, very simple function -- maybe just easier to use np.loadtxt
  # directly
  waveNum, rad = np.loadtxt(inFile, unpack=True, skiprows=nSkip)
  outDict = {'wavenumber': waveNum , 'brightness_temperature': rad, \
    'units': '(W cm-2 sr-1)/cm-1'}

  return outDict
# end readTAPE27

def readBinary(inFile, double=True):
  """
  Read LBLRTM binary file (these are special unformatted binary files,
  written in "panel" format)

  Boo...didn't work with ASTI TAPE11...

  Input
    inFile -- string, path to binary TAPE (10, 11, 12, 13)
      output_file from LBLRTM

  Output
    outWN -- float array, wavenumbers spanning spectrum
    param -- float array of ODs, radiances, fluxes, transmittances, 
      or whatever other paramter is extracted from inFile

  Keywords
    double -- boolean, is inFile in double precision? defaults to yes
  """

  import FortranFile
  from lblTools import readTape12

  outWN, param = readTape12(inFile, double=double)

  return np.array(outWN), np.array(param)
# end readBinary()

def tempIDL(inFile, fType=0, double=True):
  """
  Read in binary TAPE files (inFile), save data to IDL save files, 
  and then read them with Python for plotting. this is pretty time-
  consuming

  This is a temporary function until I figure out how to read in 
  FORTRAN binary files in Python (looks like it can be done -- 
  https://stackoverflow.com/questions/37534220/python-read-fortran-binary-file)

  Input
    inFile -- string, path to binary TAPE (10, 11, 12, 13)
      output_file from LBLRTM

  Output

  Keywords
    fType -- int, file type (radiance, transmittance, etc.; see doc in
      /project/rc/rc2/mshep/idl/patbrown/read_lbl_file_dbl.pro)
    double -- boolean, is inFile in double precision? defaults to yes
  
  """

  from scipy.io.idl import readsav
  import utils

  # write_save_file.pro is in this Git repo:
  # https://lex-gitlab.aer.com/rpernak/common_modules
  # and is considered, along with the RC_utils.py and utils.py 
  # modules, part of the RC common library
  proFile = 'write_save_file.pro'
  if not os.path.exists(proFile):
    os.symlink('externals/common/%s' % proFile, proFile)

  if double:
    proCall = \
      "write_save_file, \'%s\', file_type=%d, /dbl" % (inFile, fType)
  else:
    proCall = \
      "write_save_file, \'%s\', file_type=%d" % (inFile, fType)
  # endif double

  idlCmd = 'idl -e "%s"' % proCall
  sOut, sErr = utils.spawn(idlCmd)

  # write_save_file.pro always writes a LBLRMT_output.sav file
  # and contains the wavenum and spectrum arrays
  tempSav = 'LBLRTM_output.sav'
  idlDat = readsav(tempSav)
  waveNum, param = idlDat['wavenum'], idlDat['spectrum']
  os.remove(tempSav)
  os.remove(proFile)

  return {'wavenumber': waveNum, 'spectrum': param}
# end tempIDL()

def readTAPE11(inFile, param='Absorption'):
  """
  Read in LBLRTM binary absorption/transmittance/reflectance file

  Call
    outDict = readTAPE11(inFile)

  Input
    inFile -- string, path to TAPE11 file

  Output
    outDict -- dictionary with wavenumber, param_vals, and param_name
      key/value pairs

  Keywords
    param -- string, should be 'Absorption', 'Reflectance', or 
      'Transmittance' (case-insensitive)
  """

  goodParams = ['Absorption', 'Reflectance', 'Transmittance']
  if param.capitalize() not in goodParams:
    sys.exit('Please specify TAPE11 param (%s), returning' % goodParams)

  outDict = {'param_name': param}

  return outDict
# end readTAPE11()

def readTAPE12(inFile, param='OD', panelsize=24):
  """
  Read in LBLRTM binary monochromatic OD/radiance file

  Call
    outDict = readTAPE12(inFile)

  Input
    inFile -- string, path to TAPE12 file

  Output
    outDict -- dictionary with wavenumber, param_vals, and param_name
      key/value pairs

  Keywords
    param -- string, should be 'OD' or 'Radiance' (case-sensitive)
  """

  goodParams = ['OD', 'Radiance']
  if param not in goodParams:
    sys.exit('Please specify TAPE12 param (%s), returning' % goodParams)

  outDict = {'param_name': param}

  # read in binary spectra
  with open(inFile, "rb") as f:
    byte = f.read(panelsize)
    while byte:
      byte = f.read(panelsize)
      print byte
    # endwhile
  # endwith

  return outDict
# end readTAPE12()

def radsumRead(inFile):
  """
  Read a single RADSUM output file and return data for a given level 
  as a dictionary to be used in radsumPlot()

  Call
    outDict = radsumRead(inFile)

  Input
    inFile -- string, path to RADSUM output file

  Output
    outDict -- dictionary with the following keys 
      (with float list values):

      up_flux: upwelling flux (W/m2) as a function of wavenumber 
        and level (nLevel x nWavenumber array)
      down_flux: downwelling flux (W/m2) as a function of wavenumber 
        and level (nLevel x nWavenumber array)
      net_flux: net flux (W/m2) as a function of wavenumber and level
        (nLevel x nWavenumber array)
      heat_rate: heating rate (K/day) as a function of wavenumber 
        and level (nLevel x nWavenumber array)
      wavenumber: spectral points (cm-1) vector (1 x nWavenumber)

  Keywords
    None
  """

  inDat = open(inFile).read().splitlines()

  # initialize lists that eventually become
  waveNum1, waveNum2 = [], []
  pLevAll, upFluxAll, dnFluxAll, netFluxAll, heatRateAll = \
    ([] for x in range(5))

  for line in inDat:
    split = line.split()

    # for each band, deduce if we are processing the output or header
    try:
      # if this works, proceed to parsing the rest of the line
      iLev = int(split[0])
    except:
      # header processing -- only wanna extract wavenumber
      if len(split) > 0:
        if split[0] == 'WAVENUMBER':
          # every WAVENUMBER string occurrence implies the start of 
          # a new block of RADSUM output, which we break up into 
          # 1 x nLev vectors for each parameter, then append to 
          # *All lists that eventually become
          # (nWavenumber x nLev) arrays of output
          if 'pLev' in locals():
            # are we past the first block of output (so pLev exists)?
            # otherwise this is unnecessary
            pLevAll.append(pLev)
            upFluxAll.append(upFlux)
            dnFluxAll.append(dnFlux)
            netFluxAll.append(netFlux)
            heatRateAll.append(heatRate)

            pLev, upFlux, dnFlux, netFlux, heatRate = \
              ([] for x in range(5))
          else:
            pLev, upFlux, dnFlux, netFlux, heatRate = \
              ([] for x in range(5))
          # endif pLev len

          # there is no waveNum "array", just an nLev-element vector
          waveNum1.append(float(split[2]))
          waveNum2.append(float(split[4]))
        # endif WAVENUMBER
      # endif split len

      continue
    # end exception

    # RADSUM output processing
    pLev.append(float(split[1]))
    upFlux.append(float(split[2]))
    dnFlux.append(float(split[3]))
    netFlux.append(float(split[4]))
    heatRate.append(float(split[5]))
  # end loop over lines

  # add last output block to arrays
  pLevAll.append(pLev)
  upFluxAll.append(upFlux)
  dnFluxAll.append(dnFlux)
  netFluxAll.append(netFlux)
  heatRateAll.append(heatRate)

  outDict = {'wavenumber1': np.array(waveNum1), \
    'wavenumber2': np.array(waveNum2), \
    'level_pressure': np.array(pLevAll), \
    'up_flux': np.array(upFluxAll), \
    'down_flux': np.array(dnFluxAll), \
    'net_flux': np.array(netFluxAll), \
    'heat_rate': np.array(heatRateAll)}

  return outDict
# end radsumRead()

def rad2BT(inWN, inRad):
  """
  Radiance to Brightness Temperature conversion courtesy of 
  https://ncc.nesdis.noaa.gov/data/planck.html

  Input
    inWN -- float array, wavenumbers of spectrum (cm-1)
    inRad -- float array, associated radiance at each wavenumber
      (standard RU used in LBLRTM: W cm-2 sr-1 / cm-1)
  Output
    outBT -- float array, brightness temperatures (K)
  """

  # convert from standard RU to mW m-2 sr-1 / cm-1
  inRad /= 1e-7

  num = 1.4387752 * inWN
  denom = np.log( (1.191042e-5 * inWN**3 / inRad) + 1 )
  outBT = num/denom

  return outBT
# end rad2BT()

def colAmt2PWV(amount):
  """
  Convert accumulated molecular amounts for total path (mol/cm2) in 
  TAPE6 to precipitable water vapor (PWV, cm)

  The conversion formula was gathered from an email with Vivienne 
  Payne ("update to conversion factor for H2O") to the AER RC email
  group on 24-Jul-2006:

  pwv (cm) = \
    (column amount)(1/avogadro)(gram molec wt h2o)(1/sp density h2o)
           = [column amnt  (molec/cm^2) ] x  2.99150e-23 (cm^3/molec)

  Call
    pwv = colAmt2PWV(amount)

  Input
    amount -- float array, column amounts for a given molecule 
      (mol/cm2)

  Output
    pwv -- float array, corresponding precipitable water vapor 
      (cm)
  """

  return amount * 2.99150e-23
# end colAmt2PWV()

def wvAmtTAPE6(inTAPE6):
  """
  Extract water vapor accumulated (over entire column) amount from 
  inTAPE6
  """

  dat = open(inTAPE6).read().splitlines()
  for iLine, line in enumerate(dat):
    if 'ACCUMULATED' in line:
      wvLayAmt = float(dat[iLine+1][57:70])
    else:
      continue
    # endif ACCUMULATED

    # break out of the loop as soon after first ACCUMULATED line
    # (otherwise we will grab other molecule densities)
    break
  # end 

  return wvLayAmt
# end wvAmtTAPE6


