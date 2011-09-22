#!/bin/sh
# ems_MAKE.sh   (shamelessly copied from spm_MAKE.sh)
# 
# Compiles the C-routines ems_getMRFparams.c and ems_getMRFprior.c 
# 
# Make the same mex-related changes here as you did in spm_MAKE.sh; 
# then run "./ems_MAKE.sh PLATFORM", where PLATFORM is the same as 
# you used for running spm_MAKE.sh
#
# ------------------------------------------------------------------------
# ems_MAKE.sh         Koen Van Leemput - August 17, 2001
# 


if [ $# = 0 ]; then
	arch="sun";
else
	arch=$1;
fi


case $arch in
    sun)
	cmex5="mex     COPTIMFLAGS=-xO5";;
    windows)
	deff=-DSPM_WIN32
	cmex5="mex.bat $deff ";;
    gcc)
	cmex5="mex     COPTIMFLAGS=-O2 -f gccopts.sh";;
    sgi)
	cmex5="mex";;
    sgi64)
	cmex5="mex";;
    hpux)
	cmex5="mex     COPTIMFLAGS=-O";;
    *)
	echo "Sorry, not set up for architecture $arch"
	exit;;
esac


echo "Compiling mex files..."
$cmex5 ems_getMRFparams.c 
$cmex5 ems_getMRFprior.c 


echo "Done."
