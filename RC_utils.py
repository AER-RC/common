import os, sys
import numpy as np

class constants():
  def __init__(self, cgs=False):
    """
    Just a bunch of physical constants that we'll use at some point
    """

    # all attributes in mks units
    self.kB = 1.38064852e-23 # Boltzman's
    self.g = 9.80665 # acceleration due to gravity
    self.R = 2.87058e2 # dry-air gas constant
    self.nA = 6.022140857e23 # Avogadro's number
    self.c = 2.99792458e8 # speed of light, vaccum
    self.SB = 5.67051e-8 # Stefan-Boltzman (Wm-2K-4)
    self.sPerDay = 60 * 60 * 24 # seconds per day

    # wet and dry masses are from a KCP script check_wtot2.pro
    # they were used for RFMIP layer density calculations since we 
    # wanted to use hydrostatics and not the ideal gas law
    self.mH2O = 1.8016e-2 
    self.mDry = 2.8964e-2

    if cgs: self.mks2cgs()

  # end constructor

  def mks2cgs(self):
    self.kB = 1.38064852e-16 # Boltzman's
    self.g = 9.80665e2 # acceleration due to gravity
    self.R = 2.87058e6 # dry-air gas constant
    self.c = 2.99792458e10 # speed of light, vaccum
    self.mH2O = 1.8016e1
    self.mDry = 2.8964e1
    self.SB = 5.67051e-5 # Stefan-Boltzman

    return self
  # end mks2cgs
# end constants

def readTAPE7(inFile, xsTAPE7=False):
  """
  Read in a single TAPE7 (LBLRTM profile/layer amounts as calculated 
  by LBLATM subroutine) and return parameters in dictionary

  Call
    outDict = readTAPE7(inFile)

  Input
    inFile -- string, path to TAPE7 file

  Output
    outDict -- dictionary with the following keys:
      format: int, pressure format specification
      n_layers: int, number of layers
      n_molecules: int, number of molecules specified in profile
      scale_factor: float array, secant scaling factor (nLayers)
      end_alt: float array, instrument altitude
      obs_alt: float array, observer altitude
      SZA: float, solar zenith angle

      p_lay: float array, average pressure of given layer (nLayers)
      T_lay: float array, average temperature of given layer (nLayers)
      type_lay: int array, path type
      path_lay: int array, direction of path (nLayers)

      p_lev: float array, pressure at layer boundaries (nLevels)
      alt_lev: float array, altitudes at layer bounds (nLevels)
      T_lev: float array, temperature at layer bounds (nLevels)

      scale_factor_lay: float array, secant scaling factor for each
        layer (nLayers)
      vmr: float array, mixing ratio for each gas at each layer 
        (nMol-1 x nLayer; the -1 is for the broadening density)

  Keywords
    xsTAPE7 -- boolean, XS TAPE7 files are a bit different because t
      there are effectively 2 sections (HITRAN molecule section, XS 
      molecule section). this keyword forces outDict to contain two 
      additional keys:

        xs_names -- string array, names of all selected XS molecules
        xs_den -- float array, densities for all selected XS molecules
          (nXS x nLayer dimensions)
  """

  def stringSlice(inStr, idxArr):
    """
    Return substring of inStr that spans indices from idxArr
    """

    outStr = inStr[idxArr.min():idxArr.max()+1]

    # if it's an empty string, replace with zero
    outStr = 0 if len(outStr.strip()) == 0 else outStr

    return outStr
  # end stringSlice()

  # skip header of TAPE7, keep everything else
  datT7 = open(inFile).read().splitlines()[1:]
  record21 = datT7[0]
  profile = datT7[1:]

  iForm = int(stringSlice(record21, np.array([1])))
  nLay = int(stringSlice(record21, np.array([2, 5])))
  nMol = int(stringSlice(record21, np.array([5, 10])))
  secnto = float(stringSlice(record21, np.array([10, 19])))
  h1 = float(stringSlice(record21, np.array([40, 48])))
  h2 = float(stringSlice(record21, np.array([52, 60])))
  sza = float(stringSlice(record21, np.array([65, 73])))

  # the number of "layer lines" is dependent on the number of 
  # molecules. the convention is dictated by Record 2.1.2 in the 
  # LBLRTM instructions HTML file (8 molecules per line)
  # nLayLines = P/T line + Mixing Ratios lines (records 2.1.1 + 2.1.2)
  if nMol <= 7:
    nLayLines = 1
  else:
    # +1 for the broadener ("molecule" 8)
    nLayLines = np.ceil((nMol+1)/8.0)
  # endif nMol

  # for the P/H/T line
  nLayLines += 1 

  # how the data are read (i.e., array slicing) depends on iForm
  # this is record 2.1.1
  if iForm == 0:
    ipLay = np.array([0, 10])
    itLay = np.array([10, 21])
    iSecant = np.array([21, 30])
    iType = np.array([30, 33])
    iPath = np.array([33, 35])
    iAlt1 = np.array([36, 43])
    ipLev1 = np.array([43, 51])
    itLev1 = np.array([51, 58])
    iAlt2 = np.array([58, 65])
    ipLev2 = np.array([65, 73])
    itLev2 = np.array([73, 80])
  else:
    ipLay = np.array([0, 15])
    itLay = np.array([15, 25])
    iSecant = np.array([25, 35])
    iType = np.array([35, 38])
    iPath = np.array([38, 40])
    iAlt1 = np.array([41, 48])
    ipLev1 = np.array([48, 56])
    itLev1 = np.array([56, 63])
    iAlt2 = np.array([63, 70])
    ipLev2 = np.array([70, 78])
    itLev2 = np.array([78, 85])
  # endif iForm

  doXS = False

  # assemble lists for each profile parameter
  pLay, tLay, scaleLay, typeLay, pathLay, altLev, pLev, tLev = \
    ([] for i in range(8))

  # molecule volume mixing ratios will first be a list of lists
  # "sub" VMR is a subset of all VMRs for a given layer
  vmr, subVMR = [], []
  for iLine, line in enumerate(profile):
    if 'CROSS-SECTIONS' in line:
      break
      # only need the broadener for now, not the XS densities
      # start of cross section part of the profile
      xsLine = int(iLine)
      doXS = True
      nXS = int(line.split()[0])
      nLayLines = np.ceil((nXS+1)/8.0) + 1

      # reset all of the arrays; these guys should not change WRT
      # the HITRAN molecule section
      pLay, tLay, scaleLay, typeLay, pathLay, altLev, pLev, tLev = \
        ([] for i in range(8))
      vmr, subVMR = [], []

      # now proceed just like HITRAN section
      continue
    # endif XS

    if doXS:
      # the line after the XS header is also different from the HITRAN
      # section and contains the XS species names

      if iLine == xsLine+1:
        xsNames = line.split()
        if 'OTHER' in xsNames: xsNames.remove('OTHER')
        continue
      else:
        # treat the XS bit just like the HITRAN profile, which means 
        # we have to reset the iLine counter (while also taking into
        # account the TAPE7 XS header)
        iLine -= (xsLine + 2)
      # endif iLine
    # endif XS

    if iLine % nLayLines == 0:
      # keep the HITRAN section data instead of the XS
      if doXS: continue

      # layer P/T/Z info
      pLay.append(stringSlice(line, ipLay))
      tLay.append(stringSlice(line, itLay))
      scaleLay.append(stringSlice(line, iSecant))
      typeLay.append(stringSlice(line, iType))
      pathLay.append(stringSlice(line, iPath))

      if iLine == 0:
        # the first layer has the info for the first 2 boundaries
        altLev.append(stringSlice(line, iAlt1))
        pLev.append(stringSlice(line, ipLev1))
        tLev.append(stringSlice(line, itLev1))
      # endif iLine

      altLev.append(stringSlice(line, iAlt2))
      pLev.append(stringSlice(line, ipLev2))
      tLev.append(stringSlice(line, itLev2))

      # reset this guy every layer
      subVMR = []
    else:
      # layer molecule amounts
      subVMR += line.split()

      # are we on the last line of the layer?
      if iLine % nLayLines == nLayLines-1: vmr.append(subVMR)
      # endif nLayLines
    # end modulo 0
  # end layer loop

  # convert lists to arrays
  pLay, tLay, scaleLay, typeLay, pathLay, altLev, pLev, tLev, vmr = \
    np.array(pLay), np.array(tLay), np.array(scaleLay), \
    np.array(typeLay), np.array(pathLay), np.array(altLev), \
    np.array(pLev), np.array(tLev), np.array(vmr)

  outDict = {'n_layers': nLay, 'n_molecules': nMol, 'format': iForm, \
    'scale_factor': secnto, 'obs_alt': h1, 'end_alt': h2, 'SZA': sza}

  # extract broadening density from molecule VMR and transpose VMR
  # so it is nMol x nLay
  iBroad = 7
  broadener = vmr[:, iBroad]
  iVMR = np.delete(np.arange(nMol+1), iBroad)
  vmr = vmr[:, iVMR].T

  # make a list of all lists, loop through it, and convert all lists
  # to arrays of the proper type and stuff them into outDict
  dictKeys = ['p_lay', 'T_lay', 'scale_factor_lay', 'type_lay', \
    'path_lay', 'alt_lev', 'p_lev', 'T_lev', 'vmr', 'broadener']
  tempList = [pLay, tLay, scaleLay, typeLay, pathLay, altLev, \
    pLev, tLev, vmr, broadener]
  for iKey, temp in enumerate(tempList):
    outDict[dictKeys[iKey]] = temp.astype(float)

  return outDict
# end readTAPE7

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

def readBinary(inFile, double=False):
  """
  Read LBLRTM binary file (these are special unformatted binary files,
  written in "panel" format; and judging by 
  /project/rc/rc2/mshep/idl/patbrown/read_lbl_file.pro, these files 
  are unformatted sequential files -- see 
  https://stackoverflow.com/questions/23377274/how-to-read-fortran-77-unformatted-binary-file-into-python)
  and http://www.harrisgeospatial.com/docs/Reading_and_Writing_FORT.html

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

  import array
  from scipy.io import FortranFile

  format = np.float64 if double else np.float32
  nPoints = 2 if double else 1

  with FortranFile(inFile, "r") as inFP:
    #header = np.array(['' for i in range(80)])
    header = inFP.read_record(dtype='80c')
    print(len(header))
  # endwith
  sys.exit()

  header = ''
  iRec = np.zeros(1, dtype=np.uint32)
  with open(inFile, "rb") as inFP:
    iRec = np.fromfile(inFP, dtype='uint32', count=1)
    print(iRec)
    header = np.fromfile(inFP, dtype='c', count=80)
    print(''.join(header))
    iRec = np.fromfile(inFP, dtype='uint32', count=1)
    print(iRec)
    print

    iRec = np.fromfile(inFP, dtype='uint32', count=1)
    print(iRec)
    pv1, pv2 = np.fromfile(inFP, dtype=np.float64, count=2)
    print(pv1, pv2)
    pdv = np.fromfile(inFP, dtype=format, count=1)
    print(pdv)
    iRec = np.fromfile(inFP, dtype='uint32', count=1)
    print(iRec)
    print

    iRec = np.fromfile(inFP, dtype='uint32', count=1)
    print(iRec)
    np_arr = np.fromfile(inFP, dtype=np.long, count=nPoints)
    print(np_arr)
    iRec = np.fromfile(inFP, dtype='uint32', count=1)
    print(iRec)
    print

    sys.exit()
    od = np.fromfile(inFP, dtype=format, count=2400)
    print(len(header))
    print(od[:100])

    """
    panelHeader = array.array(format)
    panelHeader.fromfile(inFP, 3)
    print(panelHeader)

    inFP.read_record(dtype='d3')

    panelHeader = array.array(format)
    panelHeader.fromfile(f, 3)
    print(panelHeader)

    panelHeader = array.array('l')
    panelHeader.fromfile(f, np)
    print(panelHeader)
    panelHeader = array.array('d')
    panelHeader.fromfile(f, 178)
    wnDat = array.array('d')
    wnDat.fromfile(f, 3)
    wnDat = np.array(wnDat)
    numFreq = array.array('l')
    numFreq.fromfile(f, 2)
    nFreq = np.array(numFreq)[0]

    # now read in the ABSCO array, which is dependent on nFreq
    abscoArr = array.array('d')
    abscoArr.fromfile(f, nFreq)
    abscoArr = np.array(abscoArr)

    # don't need any of the rest of the "panel header" garbage
    # but we do need to make a wavenumber array associated w/ 
    # absorption coefficients
    waveNum = np.arange(wnDat[0], wnDat[1]+wnDat[2], wnDat[2])

    # save the spectrum
    abscoList.append(abscoArr)
    wnList.append(waveNum)
    """
  # endwith

  return #np.array(outWN), np.array(param)
# end readBinary()

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
      level_pressure: pressure at layer boundaries 
        (nLevel-element array)

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
    # sometimes radsum has bad pressures because of string formatting
    try:
      pLev.append(float(split[1]))
    except:
      pLev.append(np.nan)

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

def fluxToHR(flux):
  """
  Flux-to-heating rate calculation using the Stefan Boltzman law

  Input
    flux -- float array, fluxes in Wm-2
    ASSUMED TO BE NLAY X NPROFILE X NBAND! or at least NLAY is first
    dimension

  Output
    hr -- float array, corresponding heating rates in T day-1
  """

  conObj = constants()

  return (np.diff(flux, axis=0) / conObj.SB)**(1/4.) / conObj.sPerDay
# end fluxToHR()

