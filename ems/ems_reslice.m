function ems_reslice(targetName, objectNames)
%
% FORMAT ems_reslice(targetName, objectNames)
%
% Reslices objectNames to the image grid of targetName, using
% trilinear interpolation.
%
% The resliced images are written to the same directory as the
% original images, and get prefix 'r'.
%
% ------------------------------------------------------------------------
% ems_reslice.m    Koen Van Leemput - August 17, 2001

if (nargin==0)
  spm_figure('Clear','Interactive');
  targetName = spm_get(1,'.img','Select reference image');
  objectNames = spm_get(Inf, '.img', 'Select images to reslice');
  ems_reslice(targetName, objectNames);
  spm_figure('Clear','Interactive');
  return;
end

if (nargin~=2)
  error('Two input arguments required');
end

spm_reslice(str2mat(targetName, objectNames), ...
    struct('mask',0,'mean',0,'hold',1,'which',1));



% SPM completely messes up the voxel size in the header files of the
% resampled images for some wonderful reason, so let's correct for
% that

[refDIM refVOX refSCALE refTYPE refOFFSET refORIGIN refDESCRIP] = ...
    spm_hread(targetName);

for i=1:size(objectNames,1)
  thisObjectName = deblank(objectNames(i,:));
  [pth,nm,xt,vr] = fileparts(deblank(thisObjectName));
  thisResampledName = fullfile(pth, ['r' nm xt vr]);
  [DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = ...
      spm_hread(thisResampledName);
  spm_hwrite(thisResampledName, ...
	     DIM,refVOX,SCALE,TYPE,OFFSET,refORIGIN,DESCRIP);
end

