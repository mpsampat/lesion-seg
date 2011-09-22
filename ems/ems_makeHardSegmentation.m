function ems_makeHardSegmentation(classificationNames)
% FORMAT ems_hardsegment(classificationNames)
%
% Makes a hard segmentation from probability maps by assigning every
% voxel completely to the class it most probably belongs to.
%
% - classificationNames: matrix of input image filenames. These filenames 
% *must* have suffix "_segX", as generated by EMS_SEGMENT.
%
% The new images will get suffix '_hardX' replacing the suffix "_segX"
% in the probability maps.
%
% ------------------------------------------------------------------------
% ems_makeHardSegmentation.m    Koen Van Leemput - August 17, 2001


% Some user interface
if (nargin==0)
  spm_figure('Clear','Interactive');
  classificationNames = spm_get(Inf,'.img','Images to work on');
  set(spm_figure('FindWin','Interactive'), 'Name', ...
      'making hard segmentation', 'Pointer','Watch');
  ems_makeHardSegmentation(classificationNames);
  spm_figure('Clear','Interactive');
  return;
end


% Check correctness of input 
nrOfClasses = size(classificationNames ,1);

[refDIM refVOX refSCALE refTYPE refOFFSET refORIGIN] = ...
    spm_hread(deblank(classificationNames(1,:)));
refImage2world = spm_get_space(classificationNames(1,:));
hardClassificationNames = [];

for class=1:nrOfClasses
  thisClassName  = deblank(classificationNames(class,:));
  thisClassNameDir = spm_str_manip(thisClassName,'Hv');
  thisClassNameFil = [spm_str_manip(thisClassName,'stv'),'.img'];

  segInd = max(findstr(thisClassNameFil, '_seg'));
  if ~isempty(segInd)
    thisHardClassName = [thisClassNameDir '/' ...
          thisClassNameFil(1:segInd-1) '_hard' num2str(class) '.img'];
    hardClassificationNames = str2mat(hardClassificationNames, ...
        thisHardClassName);
  else
    error('All input files must have extension ''_segX''')
  end
  [DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(thisClassName);
  if (~all(DIM==refDIM) | ~all(VOX==refVOX))
    error('All input files must have the same spatial extend and resolution')
  elseif ((round(1/SCALE)~=255) | (TYPE~=2))
    error('All input files must be uint8, 255 meaning 1')
  end
  image2world = spm_get_space(thisClassName);
  if ~all(image2world(:)==refImage2world(:))
    error('All input files must have same image-to-world transformation');
  end
end

hardClassificationNames = hardClassificationNames(2:end,:);


% Read classification
classification(refDIM(1), refDIM(2), refDIM(3), nrOfClasses) = uint8(0);
for class=1:nrOfClasses
  thisClassName  = deblank(classificationNames(class,:));
  fid = fopen(thisClassName, 'r');
  if (fid==-1)
    error(['Could not read ' thisClassName])
  end
  disp(['Reading ' thisClassName])
  classification(:,:,:,class) = reshape(fread(fid, inf, 'uint8'), refDIM(1:3));
  fclose(fid);
end


% Convert classification into a hard one
fprintf('Converting classification into a hard one      ');
planeClassification(DIM(1), DIM(2), nrOfClasses) = uint8(0);
for plane=1:refDIM(3)
  for class=1:nrOfClasses
    planeClassification(:,:,class) = classification(:,:,plane,class);
  end
  
  [maxis maxind] = max(planeClassification, [], 3);
  for class=1:nrOfClasses
    classification(:,:,plane,class) = (maxind == class) ...
        .* double((maxis~=0)) * 255;
  end
  fprintf('\b\b\b\b\b');
  fprintf('%-3d %%', round(100*plane/refDIM(3)));
end
fprintf('\n');
  

% Write classification
for class=1:nrOfClasses
  thisHardClassName  = deblank(hardClassificationNames(class,:));
  fid = fopen(thisHardClassName, 'w');
  if (fid==-1)
    error(['Could not write ' thisHardClassName])
  end
  disp(['Writing ' thisHardClassName]);
  fwrite(fid, classification(:,:,:,class), 'uint8');
  fclose(fid);
  DESCRIP = 'binary segmentation';
  OFFSET = 0;
  spm_hwrite(thisHardClassName, refDIM, refVOX, refSCALE, refTYPE, ...
      OFFSET, refORIGIN, DESCRIP);
  spm_get_space(thisHardClassName, refImage2world);
end

return

