function [] = ems_defaults()
%
% Indicates where to find SPM.
% 
% This file is intended to be customized for the site.
%
% Care must be taken when modifying this file.
%
% ------------------------------------------------------------------------
% ems_defaults.m    Koen Van Leemput - August 23, 2001

global SPM_DIR MIRIT_CMD 

% Specify here where SPM99 can be found
SPM_DIR = '/home/vanleemp/matlab/spm99/';

% Multi-modal affine registration based on maximization of Mutual
% Information is implemented in ems_mireg.m. However, if you have
% access to MIRIT, the original (and faster) program developed by
% Frederik Maes (frederik.maes@uz.kuleuven.ac.be), it will be used
% instead of ems_mireg.m if you specify here where it is
% located. Otherwise, simply leave this field unchanged.
MIRIT_CMD = '';
%MIRIT_CMD = '/home/vanleemp/MIRIT/Mirit.sh';









