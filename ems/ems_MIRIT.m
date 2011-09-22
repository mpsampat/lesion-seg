function ems_MIRIT(objectName, targetName, others, free, interpol, ...
    startingParams)
% 
% FORMAT ems_MIRIT(objectName, targetName, others, free, ...
%       interpol, startingParams)
%
% Affine registration of multi-modal images with MIRIT. 
%
% MIRIT is a registration program written by Frederik Maes
% (frederik.maes@uz.kuleuven.ac.be) that searches for the affine
% transformation that maximizes the mutual information between two
% images:
%
%   F. Maes, A. Collignon, D. Vandermeulen, G. Marchal, P. Suetens,
%   Multimodality image registration by maximization of mutual
%   information , IEEE transactions on Medical Imaging, vol. 16, no. 2,
%   pp. 187-198, April 1997 (in December 2000 recognized by the IEEE as
%   the journal's most frequently cited paper published in 1997)
%
% and
%
%   F. Maes, D. Vandermeulen, P. Suetens, Comparative evaluation of
%   multiresolution optimization strategies for multimodality image
%   registration by maximization of mutual information , Medical image
%   analysis, vol. 3, no. 4, pp. 373-386, December 1999
%
%
% - objectName: object image (the image to reposition).
% - targetName: target image (the image to act as template).
% - others    : other images that will undergo the same transformation
%               as objectName  
% - free      : one for free parameters, zero otherwise [1x12]
% - interpol  : registration interpolation [1=Nearest neigbour,
%               2=Trilinear, 3=Partial volume] 
% - startingParams: starting parameters
%
% ------------------------------------------------------------------------
% ems_MIRIT.m    Koen Van Leemput - August 17, 2001


% Check input arguments
if (nargin==5)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
elseif (nargin==4)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  interpol = 2;
elseif (nargin==3)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  interpol = 2;
  free = ones(1,12);
elseif (nargin==2)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  interpol = 2;
  free = ones(1,12);
  others = '';
elseif (nargin==1)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  interpol = 2;
  free = ones(1,12);
  others = '';

  global SWD
  targetName = fullfile(SWD, 'templates/T1.img');
elseif (nargin==0)
  uiForMIRIT
  return
end

if ((size(objectName,1)~=1) | (size(targetName,1)~=1))
  error('objectName and targetName must be single images');
end


% Call upon MIRIT to calculate world-to-world transformation that
% maximizes mutual information between target and object 
startingWorld2world = spm_matrix(startingParams);
world2world = ems_callMIRIT(objectName, targetName, ...
    startingWorld2world, free, interpol);


% Now apply this world-to-world tranformation
if isempty(others)
  Images = objectName;
else
  Images   = str2mat(objectName,others);
end
for i=1:size(Images,1)
  M = spm_get_space(deblank(Images(i,:)));
  spm_get_space(deblank(Images(i,:)), world2world*M);
end

% Show some graphical output
disp('Showing some graphical output');
spm_check_registration(str2mat(targetName, objectName));
disp('Done');
spm_figure('Clear','Interactive');









% =======================================================================
% =======================================================================
% =======================================================================
function uiForMIRIT()
SPMid = spm('FnBanner',mfilename,'2.9');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','MIRIT');

nsubjects = spm_input('number of subjects',1, 'e', 1);
if (nsubjects < 1)
  spm_figure('Clear','Interactive');
  return;
end

atlasNormalization = spm_input('Normalize to atlas?', 2, 'y/n', ...
    [1; 0]);

if atlasNormalization
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  interpol = 2;
  free = ones(1,12);

  global SWD
  atlasTemplate = fullfile(SWD, 'templates/T1.img');
  for i = 1:nsubjects
    targetName = atlasTemplate;
    objectName = [];
    while size(objectName,1)<1
      objectName = spm_get(1,'.img',...
          ['select object image for subject ' num2str(i)]);
    end
    others = spm_get(Inf,'.img',...
        ['select other images for subject ' num2str(i)]);
    
    eval(['targetName'    num2str(i) ' = targetName;']);
    eval(['objectName'    num2str(i) ' = objectName;']);
    eval(['others' num2str(i) ' = others;']);
  end
else
  for i = 1:nsubjects
    targetName = [];
    while size(targetName,1)<1
      targetName = spm_get(1,'.img',...
          ['select target image for subject ' num2str(i)]);
    end
    objectName = [];
    while size(objectName,1)<1
      objectName = spm_get(1,'.img',...
          ['select object image for subject ' num2str(i)]);
    end
    others = spm_get(Inf,'.img',...
        ['select other images for subject ' num2str(i)]);
    
    eval(['targetName'    num2str(i) ' = targetName;']);
    eval(['objectName'    num2str(i) ' = objectName;']);
    eval(['others' num2str(i) ' = others;']);
  end
  
  tmp2 = [2 3 5 6 8 9 11 12];
  tmp = length(tmp2);
  nap = spm_input('# Affine Params?',3,'m',[
    ' 2 params (X & Y Translations)|'...
        ' 3 params (Translations)|'...
        ' 5 params (Translations + Pitch & Roll)|'...
        ' 6 params (Rigid Body)|'...
        ' 8 params (Rigid Body + X & Y Zoom)|'...
        ' 9 params (Rigid Body + Zooms)|'...
        '11 params (Rigid Body,  Zooms + X & Y Affine)|'...
        '12 params (Rigid Body,  Zooms & Affine)'
    ],tmp2,tmp);
  free = [ones(1,nap) zeros(1,12-nap)];
  
  interpol = spm_input('Interpolation Method?',4,'m', ...
      ['Nearest Neighbour|Trilinear Interpolation|Partial Volume'], ...
      [1 2 3], 2);
  
  startingParams = [0 0 0 0 0 0  1 1 1 0 0 0];
end


spm('Pointer','Watch');
for i=1:nsubjects
  set(spm_figure('FindWin','Interactive'),...
      'Name',['Coregistering subject ' num2str(i)],'Pointer','Watch');
  drawnow;
  eval(['targetName    =    targetName' num2str(i) ';']);
  eval(['objectName    =    objectName' num2str(i) ';']);
  eval(['others = others' num2str(i) ';']);
  
  ems_MIRIT(objectName, targetName, others, free, interpol, startingParams)
  
end
spm_figure('Clear',spm_figure('FindWin','Interactive'));
spm('FigName','MIRIT: done',Finter,CmdLine);
spm('Pointer');
return;




% =======================================================================
% =======================================================================
% =======================================================================
function world2world = ems_callMIRIT(objectName, targetName, ...
    startingWorld2world, free, interpol)


% MIRIT only works on signed 16 bit images, so convert if necessary
[dummy1 dummy2 dummy3 TYPE dummy4 dummy5 dummy6] = spm_hread(objectName);
objectConverted = 0;
if TYPE~=4
  objectName = ems_convertToShort(objectName);
  objectConverted = 1;
end

[dummy1 dummy2 dummy3 TYPE dummy4 dummy5 dummy6] = spm_hread(targetName);
targetConverted = 0;
if TYPE~=4
  targetName = ems_convertToShort(targetName);
  targetConverted = 1;
end





% It's faster to consider the lowest-resolution image as the floating image
if (prod(spm_hread(targetName)) <= prod(spm_hread(objectName)))
  fltimg = targetName;
  refimg = objectName;
else  
  fltimg = objectName;
  refimg = targetName; 
  startingWorld2world = inv(startingWorld2world);
end
  
  
% These are temporary files
thisTime = clock;
clockString = [datestr(datenum(thisTime(1), thisTime(2), thisTime(3))) ...
      '-' sprintf('%.2d', thisTime(4)) ':' sprintf('%.2d', thisTime(5)) ...
      ':' sprintf('%f', thisTime(6))];
inputfile = ['MIRIT.input-' clockString];
imagesfile = ['MIRIT.images-' clockString];
paramsfile = ['MIRIT.params-' clockString];
outputfile = ['MIRIToutput-' clockString '.m'];


% Write temporary files
succes = writeMIRITimagesFile(imagesfile, refimg, fltimg);
if ~succes
  error(['Problems with writing ' imagesfile]);
end
succes = writeMIRITparamsFile(paramsfile, startingWorld2world, free, ...
    interpol, fltimg, refimg);
if ~succes
  error(['Problems with writing ' paramsfile]);
end
succes = writeMIRITinputFile(inputfile, imagesfile, paramsfile, ...
    outputfile);
if ~succes
  error(['Problems with writing ' inputfile]);
end

% Call MIRIT
try
  global MIRIT_CMD
  eval(sprintf('!%s %s', MIRIT_CMD, inputfile));
  
  % Evaluate MIRIT's results
  fid = fopen(outputfile,'r');
  F = fread(fid); 
  fclose(fid); 
  eval(char(F')) 


  ImgToImg = [1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1] * ...
      RefToFltImgToImgTransform * ...
      [1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];
  if strcmp(fltimg, objectName)
    ImgToImg = inv(ImgToImg);
  end
  world2world = (spm_get_space(targetName)) * ImgToImg * ...
      inv(spm_get_space(objectName));
catch
  world2world = eye(4,4);
  disp('==============OOOOOPS!!!===============')
  disp('Could not retrieve MIRIT''s results; position not updated')
end
  
% Clean up
if 1
  delete(inputfile);
  delete(imagesfile);
  delete(paramsfile);
  delete(outputfile);

  if objectConverted
    delete([spm_str_manip(objectName,'rd') '.*']);
  end
  if targetConverted
    delete([spm_str_manip(targetName,'rd') '.*']);
  end
end  

return  



% =======================================================================
% =======================================================================
% =======================================================================
function succes = writeMIRITimagesFile(imagesfile, refimg, fltimg)

gantry = 0;

imgfid = fopen(imagesfile,'w');
if (imgfid==-1)
  succes = 0;
  return
end

fprintf(imgfid, '2\n');

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(refimg);
filename = refimg;
ORIENT = [1 2 3];
headerlength = OFFSET;

fprintf(imgfid, '%d ', ORIENT);
fprintf(imgfid, '\n');
fprintf(imgfid, '%d ', DIM);
fprintf(imgfid, '\n');
fprintf(imgfid, '%f ', VOX);
fprintf(imgfid, '\n');
fprintf(imgfid, '%f\n', gantry);
fprintf(imgfid, '%s\n', filename);
fprintf(imgfid, '%d\n', headerlength);
fprintf(imgfid, '%d %d ', [0 0 0; DIM-1]);
fprintf(imgfid, '\n');
fprintf(imgfid, '1.0 1.0 1.0\n');

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(fltimg);
filename = fltimg;
ORIENT = [1 2 3];
headerlength = OFFSET;

fprintf(imgfid, '%d ', ORIENT);
fprintf(imgfid, '\n');
fprintf(imgfid, '%d ', DIM);
fprintf(imgfid, '\n');
fprintf(imgfid, '%f ', VOX);
fprintf(imgfid, '\n');
fprintf(imgfid, '%f\n', gantry);
fprintf(imgfid, '%s\n', filename);
fprintf(imgfid, '%d\n', headerlength);
fprintf(imgfid, '%d %d ', [0 0 0; DIM-1]);
fprintf(imgfid, '\n');
fprintf(imgfid, '1.0 1.0 1.0\n');

succes = fclose(imgfid) + 1;
return



% =======================================================================
% =======================================================================
% =======================================================================
function succes = writeMIRITparamsFile(paramsfile, ...
    startingWorld2world, free, interpol, fltimg, refimg)

params_order = zeros(1,12);
orderNr = 1;
% rotation around z-axis is best performed before z-translation
free = [free(1:2) free(6) free(3:5) free(7:12)]; 
for index=1:12
  if (free(index) ~= 0)
    params_order(index) = orderNr;
    orderNr = orderNr + 1;
  end
end
params_order = [params_order(1:2) params_order(4:6) params_order(3) params_order(7:12)];

% MIRIT only supports only very limited image-to-world transfo's,
% whereas SPM can handle any affine image-to-world transfo's. Let's
% compensate for this by 'fooling' MIRIT (=let MIRIT calculate a part
% of the image-to-image transfo instead of world-to-world
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(fltimg);
MF_mirit_flt = zeros(4,4);
MF_mirit_flt(1:3,1:3) = diag(VOX);
MF_mirit_flt(1:3,4) = (-VOX.*(DIM-1)/2)';
MF_mirit_flt(4,4) = 1;
MF_residu_flt = spm_get_space(fltimg) * [1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1] * inv(MF_mirit_flt);
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(refimg);
MF_mirit_ref = zeros(4,4);
MF_mirit_ref(1:3,1:3) = diag(VOX);
MF_mirit_ref(1:3,4) = (-VOX.*(DIM-1)/2)';
MF_mirit_ref(4,4) = 1;
MF_residu_ref = spm_get_space(refimg) * [1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1] * inv(MF_mirit_ref);

world2world_mirit_start = inv(inv(MF_residu_flt) * startingWorld2world * MF_residu_ref);

params = transfoMatrix2params(world2world_mirit_start); 

% rotation params come before translation params in MIRIT
tmp = params_order(1:3);
params_order(1:3) = params_order(4:6);
params_order(4:6) = tmp;


% Top of 2-resolution strategy should be a subgrid of the original 
% image grid that best approaches 4x4x4 mm^3 grid  
[dummy1 VOX dummy2 dummy3 dummy4 dummy5 dummy6] = spm_hread(fltimg);
ideal = [4 4 4]./VOX;
downFactorTry = [floor(ideal); ceil(ideal)];
errors = 1-downFactorTry./(ones(2,1)*ideal);
[dummy ind] = min(abs(errors));
downFactor = ...
    [downFactorTry(ind(1),1) downFactorTry(ind(2),2) downFactorTry(ind(3),3)];
downFactor = max(downFactor, [1 1 1]);



crit_params = [4 256 256 1 1 interpol 2 2 2 downFactor];

fid = fopen(paramsfile, 'w');
if (fid==-1)
  succes = 0;
  return
end
fprintf(fid, '%f ', params);
fprintf(fid, '\n');
fprintf(fid, '%d ', params_order);
fprintf(fid, '\n');
fprintf(fid, '10 5.0 5.0 5.0 5.0 5.0 5.0 0.0 0.0 0.0 0.0 0.0 0.0\n');
fprintf(fid, '%d ', crit_params);
fprintf(fid, '\n');
fprintf(fid, '1 10.0 1e-3 1e-2 10 1e-5 5 0\n');
succes = fclose(fid) + 1;

return



% =======================================================================
% =======================================================================
% =======================================================================
function succes = writeMIRITinputFile(inputfile, imagesfile, paramsfile, ...
    outputfile)

fid = fopen(inputfile,'w');
if (fid==-1)
  succes = 0;
  return
end
fprintf(fid, '%s\n', imagesfile);
fprintf(fid, '%s\n', paramsfile);
fprintf(fid, '%d\n', 1);     % Registration mode
fprintf(fid, '%s\n', outputfile);

succes = fclose(fid) + 1;
return


% =======================================================================
% =======================================================================
% =======================================================================
function params = transfoMatrix2params(transfoMatrix) 
% Scale and skew
F = transfoMatrix(1:3,1:3);
PP = F * F';
d = det(F);


% Scale 
sy2 = PP(2,2) - (PP(2,3) * PP(2,3) / PP(3,3));
sx2 = PP(1,1) - (PP(1,2) * PP(3,3) - PP(1,3) * PP(2,3)) * (PP(1,2) * PP(3,3) - PP(1,3) * PP(2,3)) / (PP(3,3) * PP(3,3) * sy2);
sz2 = PP(3,3) - (PP(1,3) * PP(1,3) / sx2);

sx = sqrt(sx2);
sy = sqrt(sy2);
sz = sqrt(sz2);

% Skew
gx = PP(1,3)*sy/d;
gy = (PP(1,2)*PP(3,3) - PP(1,3)*PP(2,3))/PP(3,3) * sz / d;
gz = PP(2,3)/PP(3,3) * sz2 * sx / d;

% Sign of scale and skew parameters
S1 = F(1,1) / sx - (gy/sy)*F(1,2) + (gy*gz/sz)*F(1,3);
S2 = -(gx/sx)*F(3,1)+(gx*gy/sy)*F(3,2)+(1-gx*gy*gz)/sz*F(3,3);

if (S1)
  sx = sign(S1)*sx;
  sy = sign(S1)*sy;
end
if (S2)
  sy = sign(S2)*sy;
  sz = sign(S2)*sz;
end
if (d)
  sy=sign(d)*sy;
end

gx = sign(sy)*gx;
gy = sign(sz)*gy;
gz = sign(sx)*gz;


% Rotation
G = [1 gx*gz gx; gy 1 0; 0 gz 1];
S = [sx 0 0; 0 sy 0; 0 0 sz];
R = F * inv(G * S);

ry = -asin(R(1,3));
rx = asin(R(2,3)/(cos(ry)+eps));
rz = asin(R(1,2)/(cos(ry)+eps));

rx = 180/pi*rx;
ry = 180/pi*ry;
rz = 180/pi*rz;

% Translation
t = (transfoMatrix(1:3,4))';

% Result
params = [rx ry rz t sx sy sz gx gy gz];
return
