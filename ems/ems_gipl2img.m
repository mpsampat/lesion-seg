function ems_gipl2img(giplNames)
%
% FORMAT ems_gipl2img(giplNames)
%
% Reads GIPL-format images and converts them into SPM-format.
% The orientation is set automatically to:
%     X increases from Left to Right
%     Y increases from Posterior to Anterior
%     Z increases from Inferior to Superior
% if the orientation flag of the original gipl data is set
%
% ------------------------------------------------------------------------
% ems_gipl2img.m    Koen Van Leemput - August 17, 2001



if (nargin==0)
  giplNames = spm_get(Inf,'.gipl','select files to convert');
  ems_gipl2img(giplNames);
  spm_figure('Clear',spm_figure('FindWin','Interactive'));
  return
end


for giplNr=1:size(giplNames,1)
  filename = deblank(giplNames(giplNr,:));
  
  % Read header in gipl file
  fid = fopen(filename, 'r');
  dim = fread(fid, 4, 'short');
  image_type = fread(fid, 1, 'short');
  pixdim = fread(fid, 4, 'float');
  line1 = fread(fid, 80, 'char');
  matrix = fread(fid, 20, 'float');
  flag1 = fread(fid, 1, 'char');
  flag2 = fread(fid, 1, 'char');
  minvalue = fread(fid, 1, 'double');
  maxvalue = fread(fid, 1, 'double');
  origin = fread(fid, 4, 'double');
  pixval_offset = fread(fid, 1, 'float');
  pixval_cal = fread(fid, 1, 'float');
  user_def1 = fread(fid, 1, 'float');
  user_def2 = fread(fid, 1, 'float');
  magic_number = fread(fid, 1, 'uint');
  fclose(fid);
  
  % Change suffix from '.gipl' to '.img'
  goalfilename = [spm_str_manip(filename,'rd') '.img'];
  cmd = ['mv ' filename ' ' goalfilename];
  unix(sprintf('%s', cmd));

  
  % Set image-to-world transfo, depending on spatial orientation
  % Take a look at the flag flag1
  if flag1==9
    % Axial
    Orient=1;
  elseif flag1==10  
    % Coronal
    Orient=2;
  elseif flag1==11  
    % Sagital
    Orient=3;
  else 
    Orient = spm_input(['Orientation for    ' ...
          deblank(spm_str_manip(filename,'a20')) ' ?'], ...
        1,'m',['Axial|Coronal|Sagital'],[1 2 3]);
  end
  
  if Orient==1
    % Axial
    M = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];
  elseif Orient==2
    % Coronal
    M = [-1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1];      
  elseif Orient==3
    % Sagital
    M = [0 0 1 0; -1 0 0 0; 0 -1 0 0; 0 0 0 1];     
  end    
  
  
  % Set image-to-world transform such that the spatial orientation is set
  % along SPM standards and that the image centers of the atlas and the image
  % fall together

  global SWD
  atlasTemplate = fullfile(SWD, 'templates/T1.img');
  center_template = (spm_hread(atlasTemplate)-1)/2 + 1;
  center_template = [center_template(:); 1];
  V = [pixdim(1) 0 0 0; 0 pixdim(2) 0 0; 0 0 pixdim(3) 0; 0 0 0 1];
  tmp = inv(V) * inv(M) * ...
        spm_get_space(atlasTemplate) * ...
        center_template;
  ORIGIN = (dim(1:3)-1)/2 - tmp(1:3) + 1;
  image_to_world = M * V * [1 0 0 -ORIGIN(1);
                            0 1 0 -ORIGIN(2);
                            0 0 1 -ORIGIN(3);
                            0 0 0          1];
  spm_get_space(goalfilename, image_to_world);
  
  
  
  % Create header file
  P = goalfilename;  
  DIM = dim;      
  VOX = pixdim; 
  SCALE = 1; 
  TYPE = 4;
  OFFSET = 256;
  spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET);


end  
