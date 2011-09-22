function [] = ems_getMRFprior(varargin)
% 
% FORMAT MRFprior = ems_getMRFprior(classification, playing, atlas, ...
%                   lkp, G, H)
% 
% Gets the prior probability for realizations of the Potts Markov
% random field model, based on the mean field approximation as
% described in 
% 
%   K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens: Automated
%   model-based tissue classification of MR images of the brain, IEEE
%   transactions on medical imaging, vol. 18, no. 10, pp. 897-908,
%   October 1999
% 
% - classification: a 4-D uint8 matrix of size [I J K N], describing
%   the probability (x 255) that voxel (i,j,k) in a 3-D image grid of
%   size [I J K] belongs to a class n. Note that sum(classification,4)
%   should always equal ones(I, J, K).
% - playing: a 3-D uint8 matrix of size [I J K] that indicates whether
%   voxel (i,j,k) belongs to the region of interest (if 1) or not (if 0)
% - atlas: a 4-D uint8 matrix of size [I J K M] describing the
%   probability (x 255) that voxel (i,j,k) in a 3-D image grid of size
%   [I J K] belongs to a tissue type N. Note that sum(atlas,4) should
%   always equal ones(I, J, K).
% - lkp: look-up table of size [1 N] that describes which
%   classification maps are considered to form together one single
%   tissue type in the atlas.
% - G and H: matrices of size [M M] that contain the MRF transition
%   costs, in the plane and out of the plane, respectively. (see paper)
%
% - MRFprior: calculated prior in the same format as classification
%
% ------------------------------------------------------------------------
% ems_getMRFprior.m    Koen Van Leemput - August 17, 2001



% If you get here, then you have not compiled the C-file
error('You don''t seem to have compiled ems_getMRFprior.c yet...')

