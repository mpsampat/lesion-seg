function correctedDataNames = ems_getBiasCorrectedImage(dataNames)
%
% FORMAT correctedDataNames = ems_getBiasCorrectedImage(dataNames)
%
% Corrects images for their bias field estimated by EMS_SEGMENT or
% EMS_LESIONS. The resulting image names get suffix '_corr'.
%
% ------------------------------------------------------------------------
% ems_getBiasCorrectedImage.m    Koen Van Leemput - August 17, 2001


if (nargin==0)
  SPMid = spm('FnBanner',mfilename,'2.9');

  dataNames = spm_get(Inf,'.img', 'Select original images');
  
  ems_getBiasCorrectedImage(dataNames);

  spm_figure('Clear',spm_figure('FindWin','Interactive'));

  return;
end  
  
toEvaluate = 'i1 ./ (i2 + (i2==0)) .* (i2~=0)';

correctedDataNames = [];
for i=1:size(dataNames,1)
  thisDataName = dataNames(i,:);
  thisBiasName = [spm_str_manip(thisDataName,'rd') '_bias.img'];

  if (exist(thisBiasName, 'file')==0)
    disp([thisBiasName ' does not exist: skipping bias field correction'])
  else
    thisCorrectedName = [spm_str_manip(thisDataName,'rd') '_corr.img']; 
    disp(['Correcting ' thisDataName])
    spm_imcalc_ui(str2mat(thisDataName, thisBiasName), ...
        thisCorrectedName, toEvaluate);
  end
  correctedDataNames = strvcat(correctedDataNames, thisCorrectedName);
end
