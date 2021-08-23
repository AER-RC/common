#!/bin/tcsh

# for implementation, see local version build directions in 
# https://github.com/AER-RC/LBLRTM/wiki/LBLRTM-New-Release-Procedure

date

# 1 command line (CLI) argument necessary
if ($#argv == 0) then
  echo "Usage: "$0 "replace"
  echo "   or: "$0 "update"
  exit 1
endif

# for cron jobs, we need the PATH specified such that the 
# libraries for pgf90 and ifort are viewable
# probably don't need all of these
set path = (/usr/lib64/qt-3.3/bin /usr/kerberos/sbin /usr/kerberos/bin \
    /usr/local/bin /usr/bin /bin /usr/local/grads/bin \
    /usr/local/ncarg/bin /nas/local/pgi_15.3/linux86-64/15.3/bin \
    /nas/local/pgi/linux86-64/10.5/bin /nas/local/intel_2018/bin \
    /nas/local/pgi_15.3/linux86-64/15.3/bin \
    /nas/local/pgi/linux86-64/10.5/bin)

set LD_LIBRARY_PATH = (/nas/local/intel_2018/lib/intel64_lin)

# AER local versions
setenv LOCAL /nas/project/rc/rc1
setenv LBLRTM lblrtm_local_version
setenv CNTNM cntnm_local_version
setenv MONO monortm_local_version
setenv LNFL lnfl_local_version
setenv LFILE line_parameters_local_version
setenv XS xs_files_local_version
setenv RS radsum_local_version

# for LFILE, should I use what i have or line_parameters_incoming/?

if ($1 == "replace") then 
  echo 'Replacing local version'

  # wipe out the current local version
  cd $LOCAL
  rm -rf $LBLRTM
  rm -rf $CNTNM
  rm -rf $MONO
  rm -rf $XS
  rm -rf $LNFL
  rm -rf $RS
  git clone --recursive git@github.com:AER-RC/LNFL.git
  git clone --recursive git@github.com:AER-RC/LBLRTM.git
  git clone --recursive git@github.com:AER-RC/mt-ckd.git
  git clone --recursive git@github.com:AER-RC/cross-sections.git
  git clone --recursive git@github.com:AER-RC/monoRTM.git
  git clone --recursive git@github.com:AER-RC/RADSUM.git

  # stopped version controlling Line File at v3.6 because the size 
  # of the repo was unsustainable
  #cd $LFILE
  #rm -rf AER_line_files
  #git clone git@lex-gitlab.aer.com:rpernak/AER_Line_File.git AER_line_files

  # stopped giving others write permission in $LOCAL because development
  # SHOULD BE HAPPENING IN BRANCHES AND NOT IN THE LOCAL VERSION
  # local versions exist only for exectuables that all can use on 
  # AER systems
  #cd $LOCAL

  #chmod -R g+w $LBLRTM
  #chmod -R g+w $CNTNM
  #chmod -R g+w $MONO
  #chmod -R g+w $XS
  #chmod -R g+w $RS
  #chmod -R g+w $LNFL
  #chmod -R g+w $LFILE
else if ($1 == "update") then
  echo 'Updating local version'

  # update instead of checkout
  cd $LOCAL/$LBLRTM
  echo "LBLRTM update"
  git pull
  git submodule update --init --recursive

  cd $LOCAL/$CNTNM
  echo "MT_CKD update"
  git pull
  git submodule update --init --recursive

  cd $LOCAL/$MONO
  echo "monoRTM update"
  git pull
  git submodule update --init --recursive

  cd $LOCAL/$LNFL
  echo "LNFL update"
  git pull
  git submodule update --init --recursive

  cd $LOCAL/$RS
  echo "radsum update"
  git pull
  git submodule update --init --recursive

  cd $LOCAL/$XS
  echo "XS update"
  git pull

  #cd $LOCAL/$LFILE/AER_line_files
  #git pull
else
  # invalid CLI argument
  echo 'Usage: '$0' replace'
  echo '   or: '$0' update'
  exit 1
endif # replace/update

# rebuild executables
cd $LOCAL/$CNTNM/build
rm ../cntnm_v*
gmake -f make_cntnm linuxPGIdbl
gmake -f make_cntnm linuxPGIsgl
gmake -f make_cntnm linuxINTELdbl
gmake -f make_cntnm linuxINTELsgl
rm -rf *.mod cntnm_v*

cd $LOCAL/$LBLRTM/build
rm ../lblrtm_v*
gmake -f make_lblrtm linuxPGIdbl
gmake -f make_lblrtm linuxPGIsgl
gmake -f make_lblrtm linuxINTELdbl
gmake -f make_lblrtm linuxINTELsgl
rm -rf *.mod lblrtm_v*

cd $LOCAL/$LNFL/build
rm ../lnfl_v*
#gmake -f make_lnfl linuxPGIdbl
gmake -f make_lnfl linuxPGIsgl
gmake -f make_lnfl linuxINTELsgl
rm -rf lnfl_v*
# skip monoRTM for now
#exit 0

cd $LOCAL/$MONO/src
# if we don't do these moves, the builds will be attempted with the 
# netCDF option, which we decided not to use
# "In order to generate the netcdf files containing the ODs by layer, 
# frequency and species the user sets a flag when running make."
mv -f netcdf_helper_mod.F90 netcdf_helper_mod.F90.safe
mv -f netcdf_helper_mod.f90 netcdf_helper_mod.f90.safe

cd $LOCAL/$MONO/build
# don't know what IBM machine to use, so these builds won't work
#gmake -f make_monortm aixIBMdbl
#gmake -f make_monortm aixIBMsgl
rm ../monortm_v*
gmake -f make_monortm linuxPGIsgl
gmake -f make_monortm linuxPGIdbl
gmake -f make_monortm linuxINTELsgl
gmake -f make_monortm linuxINTELdbl
rm -rf monortm_v*

cd $LOCAL/$RS/build
rm ../radsum_v*
gmake -f make_radsum linuxPGIdbl
gmake -f make_radsum linuxINTELdbl
rm -rf radsum_v*

exit 0

