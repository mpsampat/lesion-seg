function ems_compareToAtlas(dataNames)
%
% FORMAT ems_compareToAtlas(dataNames)
%
% Check spatial position of images with respect to the atlas 
%
% ------------------------------------------------------------------------
% ems_compareToAtlas.m    Koen Van Leemput - August 17, 2001


if (nargin==0)
  SPMid = spm('FnBanner',mfilename,'2.9');

  dataNames = spm_get(Inf,'.img', 'Select images to compare with atlas');
  
  ems_compareToAtlas(dataNames);

  spm_figure('Clear',spm_figure('FindWin','Interactive'));

  return;
end  




global SWD
atlasName = fullfile(SWD, 'templates/T1.img');

spm_check_registration(str2mat(atlasName, dataNames));
