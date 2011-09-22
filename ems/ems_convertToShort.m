function convertedDataNames = ems_convertToShort(dataNames)
%
% FORMAT convertedDataNames = ems_convertToShort(dataNames)
%
% Converts images into spm type 4 (signed 16 bit). 
%
% The resulting images are written with '_short' appended to the
% original file names.
%
% ------------------------------------------------------------------------
% ems_convertToShort.m     Koen Van Leemput - August 17, 2001

if (nargin==0)
  dataNames = spm_get(Inf,'.img','select files to convert');
  ems_convertToShort(dataNames);
  spm_figure('Clear',spm_figure('FindWin','Interactive'));
  return
end


convertedDataNames = [];
for dataNr=1:size(dataNames,1)
  dataName = deblank(dataNames(dataNr,:));
  outputName = [spm_str_manip(dataName,'rd') '_short.img'];
  disp(['Converting ' dataName ' to int16'])
  spm_imcalc_ui(dataName, outputName, 'i1');
  convertedDataNames = strvcat(convertedDataNames, outputName);
end  


