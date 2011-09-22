Installation
============

First download the SPM software from the SPM web site
(http://www.fil.ion.ucl.ac.uk/spm/). At the same time, get the latest
patches as well (SPM bugs & fixes area:
ftp://ftp.fil.ion.ucl.ac.uk/spm/spm99_updates). SPM installs out of
the box on some platforms, but in general you'll have to compile some
C-routines using 'spm_MAKE.sh' (you might have to edit 'spm_MAKE.sh'
to make it work on your local system).

Uncompress and untar the EMS software package. This will create a
directory 'ems' that contains the EMS functions and routines. The
first thing to do is to edit the file 'ems_defaults.m' to provide the
correct path to SPM on your system. Subsequently, run 'ems_MAKE.sh'
(if you had to edit 'spm_MAKE.sh' while installing SPM, you'll have to
make the same changes here). Start Matlab, make sure that the ems
directory is in your MATLABPATH (type 'path' at the Matlab prompt to
see your current path, and use the 'addpath' command if necessary),
and type 'ems' at the Matlab prompt. This will invoke SPM and
additionally launch the EMS menu window shown below. Installation
completed.

Getting started
===============

The EMS web site
(http://bilbo.esat.kuleuven.ac.be/web-pages/downloads/ems/) contains a
detailed description on how to use the software.

