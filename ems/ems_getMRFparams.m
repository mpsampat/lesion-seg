function [] = ems_getMRFparams(varargin)
% 
% FORMAT [G, H] = ems_getMRFparams(classification, lkp, sampleInd, ...
%                                  constraintType) 
%
% Estimates the parameters of a Potts Markov random field model
% from a classification as described in
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
% - lkp: look-up table of size [1 N] that describes which
%   classification maps are considered to form together one single
%   tissue type.
% - sampleInd: list of those voxels that are used to estimate the MRF
%   parameters. A voxel (i,j,k) is indexed as i*j*k.
% - constraintType: determines what sort of constraints are applicable
%   on the estimated MRF parameters G and H:
%          0 -> no constraints
%          1 -> adjacent classes have the same prior probability in 
%               each other's neighborhood (see paper)
%          2 -> transition matrices G and H are symmetric
%          3 -> all transition costs are the same
%
% - G and H: estimated MRF transition costs, in the plane and out
%   of the plane, respectively. (see paper)
%
% ------------------------------------------------------------------------
% ems_getMRFparams.m    Koen Van Leemput - August 17, 2001



% If you get here, then you have not compiled the C-file
error('You don''t seem to have compiled ems_getMRFparams.c yet...')
