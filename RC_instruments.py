import os, sys
import numpy as np

# should be in same directory as this module
import utils

def srfAIRS(fileSRF='/nas/ASR/RC_Release_Benchmark_Tests/AIRS/' + \
  'data/srftables_051118v4.h5'):

  """
  Read in and return the AIRS spectral response function from HDF5 
  file (amongst other parameters -- outDict output)

  Input
    None

  Output
    outDict -- dictionary with following keys:
      wavenumber -- float array, wavenumbers of spectrum (NOT the 
        line centers; rather, the wavenumber associated with each 
        point in the SRF)
      FWHM -- float array, FWHM at each wavenum
      SRF -- float array (nWN x nFGrid)
      center_line -- float array, line centers (cm-1)

  Keywords
    fileSRF -- string, full path to static spectral response function 
      (SRF) HDF5 file. this was converted from HDF4 with h4toh5 
      (https://ftp.hdfgroup.org/h4toh5/download.html)

  """

  import h5py

  utils.file_check(fileSRF)

  h5Obj = h5py.File(fileSRF, 'r')
  center = h5Obj['freq'][:]
  srf = h5Obj['srfval'][:]
  fwGrid = h5Obj['fwgrid'][:] # not entirely sure what this is...
  width = h5Obj['width'][:]
  h5Obj.close()

  # make arrays of repeated vectors such that (nWN x 1) and 
  # (1 x nFGrid) can be added and multiplied to produce a full grid 
  # of wavenumbers associated with the SRF
  nWN, nFGrid = srf.shape
  centerArr = np.tile(center, (1, nFGrid))
  widthArr = np.tile(width, (1, nFGrid))
  fwArr = np.tile(fwGrid, (nWN, 1))
  waveNum = centerArr + fwArr * widthArr

  return {'wavenumber': waveNum, 'FWHM': width, 'SRF': srf, \
    'center_line': center}
# end srfAIRS()

IDXNPZ = '/nas/ASR/RC_Release_Benchmark_Tests/AIRS/output_files/' + \
  'AIRS_LBLRTM_Overlap_Idx.npz'

def airsSRFInterpOverlap(srfCenter, srfWN, srf, modWN, \
  outNPZ=IDXNPZ):
  """
  Designed to help save time in interpolateSRF(), which loops over 
  all line centers that are in the LBLRTM spectral range and extracts
  indices of points in the line center +/- FWHM range. These should 
  stay constant with LBLRTM and Line File releases, because the AIRS 
  SRF does not change and we should be keeping the same spectral 
  ranges and resolutions with LBLRTM runs that are done in 
  benchmark_tests.py

  Input
    srfCenter -- vector array, 1 x nWavenumbers array of line centers
      for Spectral Response Function (SRF)
    srfWN -- array, nFWHM x nWavenumbers complete wavenumber array 
      for SRF (i.e., not just line centers)
    srf -- float array, corresponding spectral response function 
      values at wavenumbers in srfWN
    modWN -- vector array, 1 x nWNLBL of wavenumbers in 
      LBLRTM spectrum

  Output
    overlapIdx -- array list where each element is a list of indices 
      for a given line center (in the LBLRTM spectral range) that 
      falls in the line center +/- FWHM range. lists do not 
      necessarily contain the same amount of elements

  Keywords
    outNPZ -- string, full path to compressed NumPy file in which 
      overlapIdx is saved
  """

  from scipy.interpolate import interp1d

  # what part of AIRS spectrum (line centers) is contained by 
  # the LBL spectrum? --> indices for center overlap
  icOverlap = np.where(\
    (srfCenter >= modWN.min()) & (srfCenter <= modWN.max()))[0]

  # for each matching center line, zoom in on the center +/- FWHM 
  # region of interest (ROI)
  # pretty time consuming...
  # indices for full overlap
  ifOverlap = []; iSRF = []
  for i, iOver in enumerate(icOverlap):
    print('{:d} of {:d}'.format(i, icOverlap.size))
    wnOver = srfWN[iOver, :]

    # grab LBL spectral points in the ROI
    # this is the time hog, and neither method is faster
    #iLo = lblWN >= wnOver.min()
    #iHi = lblWN <= wnOver.max()
    #iOverLBL = iLo & iHi
    iOverLBL = \
      np.where((modWN >= wnOver[0]) & (modWN <= wnOver[-1]))[0]

    if iOverLBL.size == 0: continue

    ifOverlap.append(iOverLBL)
    wnOver = srfWN[iOver, :]
    srfOver = srf[iOver, :]

    # fit a quadratic to the AIRS spectrum, then solve equation with 
    # LBL spectral points (i.e., interpolate onto LBL grid)
    quadFunc = interp1d(wnOver, srfOver, kind='quadratic')
    interpSRF = quadFunc(modWN[iOverLBL])
    interpSRF /= sum(interpSRF)
    iSRF.append(interpSRF)
  # end iOver loop

  np.savez(\
    '%s' % outNPZ[:-4], interp_srf=iSRF, full=ifOverlap, \
      center=icOverlap)

  return icOverlap, ifOverlap, iSRF
# end airsSRFInterpOverlap()

def interpolateSRF(inTAPE12, outNPZ='temp.npz', idxNPZ=IDXNPZ):
  """
  Interpolate the AIRS spectral response function onto the grid given 
  by an LBLRTM TAPE12 spectrum. Right now this is a very time 
  consuming function because of the np.where() calls in a loop with 
  > 2000 iterations

  Call
    outWN, outRad, outBT = interpolateSRF(inTAPE12)

  Input
    inTAPE12 -- string, full path to TAPE12 from LBLRTM that contains 
      radiances for a given profile

  Output
    outWN -- float array, line centers (cm-1) of AIRS spectrum
    outRad -- float array, corresponding AIRS radiances (interpolated
      onto the LBLRTM spectral grid, [W cm-2 sr-1 / cm-1] units)
    outBT -- float array, corresponding AIRS brightness temperatures 
      ([K] units)

  Keywords
    outNPZ -- string, full path to file in which the output paramters 
      are saved for later usage (compressed NumPy (.npy) file)
    idxNPZ -- string, .npz file generated with airsSRFInterpOverlap() 
      to save time. By default, this is not provided and 
      airsSRFInterpOverlap() is run
  """

  import time

  # also should be in same directory as this module
  import RC_utils as RC

  print('Reading in AIRS SRF')
  airsDict = srfAIRS()
  airsWN = airsDict['wavenumber']
  airsCenter = airsDict['center_line'][:, 0]
  airsSRF = airsDict['SRF']

  print('Reading in AIRS LBLRTM TAPE12')
  lblWN, lblSpec = RC.readBinary(inTAPE12, double=True)

  if idxNPZ is None:
    # find line center overlap
    print('Indice file not set, running airsSRFInterpOverlap()')
    lcOver, fullOver, interpSRF = \
      airsSRFInterpOverlap(airsCenter, airsWN, airsSRF, lblWN)
  else:
    # grab line center overlap
    if not os.path.exists(idxNPZ):
      sys.exit('Could not find %s, you may have to run %s again' % \
        (idxNPZ, 'airsSRFInterpOverlap()') )

    print('Loading {}'.format(idxNPZ))
    npzDat = np.load(idxNPZ, allow_pickle=True, encoding='bytes')
    lcOver = npzDat['center']
    fullOver = npzDat['full']
    interpSRF = npzDat['interp_srf']
  # endif idxNPZ

  # for each matching center line, zoom in on the center +/- FWHM 
  # region of interest (ROI)
  print('Scaling AIRS SRF with LBLRTM radiances')
  outRad = []
  for i, iOver in enumerate(lcOver):
    print('{:d} of {:d}'.format(i+1, lcOver.size))
    iOverLBL = fullOver[i]
    tempSRF = interpSRF[i] * lblSpec[iOverLBL]
    outRad.append(sum(tempSRF))
  # end loop
  
  outWN = np.array(airsCenter)
  outRad = np.array(outRad)
  outBT = RC.rad2BT(outWN, outRad)

  # save to a compressed NumPy file for later use; savez attaches the
  # .npz extension
  np.savez(outNPZ[:-4], wavenumber=outWN, radiance=outRad, BT=outBT)

  return outWN, outRad, outBT
# end interpolateSRF()

