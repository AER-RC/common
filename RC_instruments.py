import numpy as np

# should be in same directory as this module
import utils

def srfAIRS(fileSRF='/nas/ASR/RC_Release_Benchmark_Tests/AIRS/' + \
  'data/srftables_051118v4.h5', radiance=False):

  """
  Read in the AIRS spectral response function from HDF5 file and 
  return corresponding brightness temperature (or radiance) spectrum

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
    fileSRF -- string, path to static spectral response function 
      (SRF) HDF5 file. this was converted from HDF4 with h4toh5 
      (https://ftp.hdfgroup.org/h4toh5/download.html)

    radiance -- boolean, return radiance instead of BT spectrum
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

def interpolateSRF(inTAPE12, outNPZ='temp.npz'):
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
      are saved for later usage
  """

  # also should be in same directory as this module
  import RC_instruments as RCI
  import RC_utils as RC

  print 'Reading in AIRS SRF'
  airsDict = RCI.srfAIRS()
  airsWN = airsDict['wavenumber']
  airsCenter = airsDict['center_line'][:, 0]
  airsSRF = airsDict['SRF']

  print 'Reading in AIRS LBLRTM TAPE12'
  lblWN, lblSpec = RC.readBinary(inTAPE12, double=True)

  # what part of AIRS spectrum is contained by LBL spectrum?
  iOverlap = np.where(\
    (airsCenter >= lblWN.min()) & (airsCenter <= lblWN.max()))[0]

  # for each matching center line, zoom in on the center +/- FWHM 
  # region of interest (ROI)
  # pretty time consuming...
  outRad = []
  for i, iOver in enumerate(iOverlap):
    print '%d of %d' % (i, iOverlap.size)
    wnOver = airsWN[iOver, :]
    srfOver = airsSRF[iOver, :]

    # grab LBL spectral points in the ROI
    # this is the time hog, and neither method is faster
    iLo = lblWN >= wnOver.min()
    iHi = lblWN <= wnOver.max()
    iOverLBL = iLo & iHi
    #iOverLBL = np.where((lblWN >= wnOver[0]) & (lblWN <= wnOver[-1]))[0]
    #if iOverLBL.size == 0: continue
    if lblWN[iOverLBL].size == 0: continue

    # fit a quadratic to the AIRS spectrum, then solve equation with 
    # LBL spectral points (i.e., interpolate onto LBL grid)
    quadFit = np.polyfit(wnOver, srfOver, 2)
    quadFunc = np.poly1d(quadFit)
    interpSRF = quadFunc(lblWN[iOverLBL])
    interpSRF /= sum(interpSRF)
    interpSRF *= lblSpec[iOverLBL]
    outRad.append(sum(interpSRF))
  # end iOver loop

  outWN = np.array(airsCenter)
  outRad = np.array(outRad)
  outBT = RC.rad2BT(outWN, outRad)

  # save to a compressed NumPy file for later use; savez attaches the
  # .npz extension
  np.savez(outNPZ[:-4], wavenumber=outWN, radiance=outRad, BT=outBT)

  return outWN, outRad, outBT
# end interpolateSRF()

