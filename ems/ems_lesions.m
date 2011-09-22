function [inputArguments, historyOfParams] = ems_lesions(dataNames, ...
						  use2D, maxBiasOrder, ...
						  useMRF, updateMRFparams, ...
						  G, H, mahalanobis, ...
						  lesionConstraint, ...
						  lesionHomeClass, ...
						  garbageConstraint, ...
						  garbageHomeClass)
%
% FORMAT [inputArguments, historyOfParams] = ems_lesions(dataNames, ...
%               use2D, maxBiasOrder, useMRF, updateMRFparams. G, H, ...
%               mahalanobis, lesionConstraint, lesionHomeClass, ...
%               garbageConstraint, garbageHomeClass)
%
% Automated segmentation of intensity abnormalities (such as for
% instance Multiple Sclerosis lesions) in MR images of the brain by
% detection of model outliers, as described in
%
%   K. Van Leemput, F. Maes, D. Vandermeulen, A. Colchester,
%   P. Suetens, Automated segmentation of multiple sclerosis lesions
%   by model outlier detection , IEEE transactions on medical imaging,
%   vol. 20, no. 8, pp. 677-688, August 2001
%
% The method is an extension of the model-based tissue classification
% of normal MR brain images implemented in EMS_SEGMENT.
%
%
% - dataNames: string matrix that contains in every row the filename
%   of the row-th channel of the multi-spectral MR data to be segmented.
% - use2D: if 0, use a 3D polynomial bias field model; otherwise, use
%   a 2D polynomial for each slice separately. Default is 0.
% - maxBiasOrder: order of the bias field polynomial model.  Default
%   is 4.
% - useMRF: if 0, classify each voxel only based on its intensity and
%   on the information provided by the atlas; otherwise, use an
%   additional Markov random field prior model for the
%   classification. Default is 0.
% - updateMRFparams: if 0, the Markov random field parameters G and H
%   are kept fixed at the values that were provided by the user;
%   otherwise, these parameters are automatically estimated. Default
%   is 1.
% - G and H are the Markov random field interaction matrices (4x4),
%   for the classes labeled as follows: white matter = class 1, gray
%   matter = class 2, csf = class 3, other = class 4. Default is
%   zeros(4,4).
% - mahalanobis: Mahalanobis distance at which intensities are
%   considered abnormal with respect to the normal distributions (see
%   paper)
% - lesionConstraint: string describing the intensity constraint on
%   candidate lesions. If a voxel does not satisfy this condition, it
%   will not be considered as belonging to a lesion. The format is as
%   follows: 'i1', 'i2', 'i3' etc refers to the bias-corrected
%   intensity in MR-channel 1, 2, 3 respectively. Furthermore, 'wm1',
%   'gm1', 'csf1' refers to the estimated mean intensity of white
%   matter, grey matter and csf in channel 1 ('wm2', 'gm2', 'csf2'
%   refer to the same in channel 2, etc). In the paper, where MS
%   lesions were segmented from multi-spectral data (T1-weighted,
%   T2-weighted and PD-weighted), the following expression was used:
%   '(i2>=gm2) & (i3>=gm3)', indicating that only voxels with a higher
%   intensity than the mean intensity of grey matter in the T2- and
%   PD-weighted channel were considered as candidate
%   lesions. Default is '' (no intensity constraints).
% - lesionHomeClass: string describing to which major tissue type
%   lesions are considered to belong. This has only effect when the
%   Markov random field prior is used (see paper). Choice is 'wm',
%   'gm', 'csf', corresponding to white matter, grey matter and CSF
%   respectively. Since Multiple Sclerosis lesions are typically
%   located inside white matter, this was 'wm' in the paper. Default
%   is 'wm'.
% - garbageConstraint: similar to 'lesionConstraint', but for
%   occasional other model outliers than the ones of interest
%   (somewhat unrespectfully labeled as 'garbage'). The explicit
%   detection of model outliers that are actually not of interest may
%   be important for accurate estimation of the model parameters. In
%   the paper, there were gross dark outliers in the data, so the
%   following expression was used: '(i2<gm2) & (i3<gm3)'. Default is
%   'i1>Inf' (no detection of other outliers than the ones of
%   interest). IMPORTANT: 'lesionConstraint' and 'garbageConstraint'
%   should ALWAYS be mutually exclusive, i.e. no voxel can exist that
%   satisfies both constraints simultaneously (see paper)
% - garbageHomeClass: similar to 'lesionHomeClass', but for the
%   'garbage' class. In the paper, the dark outliers were located
%   inside the CSF, so 'csf' was used. Default is 'csf'.
%
% 
% - inputArguments: contains the input arguments with which the
%   program was called
% - historyOfParams: contains the estimated model parameters at every
%   iteration.
%
% General remarks:
% - see also EMS_SEGMENT 
% - make sure you know what you are doing when using this software
%   with your own home-made intensity constraints and home classes...
% 
% ------------------------------------------------------------------------
% ems_lesions.m    Koen Van Leemput - August 17, 2001



% Get input data, check and correct if necessary
if (nargin<12)
  garbageHomeClass = 'csf';
end
if (nargin<11)
  garbageConstraint = 'i1>Inf'; 
end
if (nargin<10)
  lesionHomeClass = 'wm';
end
if (nargin<9)
  lesionConstraint = '';
end
if (nargin==7)
  mahalanobis = 3;
elseif (nargin==6)
  mahalanobis = 3;
  H = zeros(4,4);
elseif (nargin==5)
  mahalanobis = 3;
  G = zeros(4,4);
  H = zeros(4,4); 
elseif (nargin==4)
  mahalanobis = 3;
  if useMRF
    updateMRFparams = 1;
  else
    updateMRFparams = 0;
  end
  G = zeros(4,4);
  H = zeros(4,4);
elseif (nargin==3)
  mahalanobis = 3;
  useMRF = 0;
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
elseif (nargin==2)
  mahalanobis = 3;
  useMRF = 0;
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
  maxBiasOrder = 4;
elseif (nargin==1)
  mahalanobis = 3;
  useMRF = 0;
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
  maxBiasOrder = 4;
  use2D = 0;
elseif (nargin==0)
  uiForems_lesions
  return
end

if ~useMRF
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
end
checkConsistency(dataNames, use2D, maxBiasOrder, useMRF, updateMRFparams, ...
		 G, H, mahalanobis, lesionConstraint, lesionHomeClass, ...
		 garbageConstraint, garbageHomeClass)
prettyEcho(dataNames, use2D, maxBiasOrder, useMRF, updateMRFparams, ...
	   G, H, mahalanobis, lesionConstraint, lesionHomeClass, ...
	   garbageConstraint, garbageHomeClass)



% Some initialization and allocation of variables
scaleFactor=100;
lambda = exp(-.5*mahalanobis^2);
dumping = [1 1 1 0 0 0];

global SWD
atlasDirectory = fullfile(SWD, 'apriori');
atlasNames = str2mat(fullfile(atlasDirectory, 'white.img'),...
                     fullfile(atlasDirectory, 'gray.img'),...
                     fullfile(atlasDirectory, 'csf.img'));

nc = [1 1 1 3];
lkp    = []; for i=1:size(nc,2), lkp = [lkp ones(1,nc(i))*i]; end;

infoData = spm_vol(dataNames);
infoAtlas = spm_vol(atlasNames);

nrOfChannels = size(dataNames,1);
nrOfPriors = size(atlasNames,1) + 1;
nrOfClasses = size(lkp,2);

scaleData = zeros(1, nrOfChannels);
offsetData = zeros(1, nrOfChannels);
for channel=1:nrOfChannels
  scaleData(channel) = infoData(channel).pinfo(1);
  offsetData(channel) = infoData(channel).pinfo(2);
end

DIM = infoData(1).dim(1:3);
nrOf2DElements = DIM(1)*DIM(2);
nrOf3DElements = nrOf2DElements * DIM(3);

means = zeros(nrOfChannels, nrOfClasses);
covariances = zeros(nrOfChannels, nrOfChannels, nrOfClasses);
classWeights = zeros(1, nrOfClasses);
backgroundVariance = zeros(nrOfChannels,1);

inputArguments = struct(...
    'dataNames', dataNames, ...
    'maxBiasOrder', maxBiasOrder, ...
    'use2D', use2D, ...
    'useMRF', useMRF, ...
    'updateMRFparams', updateMRFparams, ...
    'G', G, ...
    'H', H, ...
    'mahalanobis', mahalanobis, ...
    'lesionConstraint', lesionConstraint, ...
    'lesionHomeClass', lesionHomeClass, ...
    'garbageConstraint', garbageConstraint, ...
    'garbageHomeClass', garbageHomeClass);



% In order to convert data into uint16-format, calculate minimum
% and maximum value first
fprintf('Calculating maximum and minimum voxel values      ');
maxVals = repmat(-Inf, [nrOfChannels 1]);
minVals = repmat(Inf, [nrOfChannels 1]);
for plane=1:DIM(3)
  for channel=1:nrOfChannels
    A = spm_matrix([0 0 -plane]);
    B = A * inv(infoData(1).mat) * infoData(channel).mat;
    img = spm_slice_vol(infoData(channel), inv(B), DIM(1:2), 1);
    maxVals(channel)  = max([max(img(:)) maxVals(channel)]);
    minVals(channel)  = min([min(img(:)) minVals(channel)]);
  end
  fprintf('\b\b\b\b\b');
  fprintf('%-3d %%', round(100*plane/DIM(3)));
end
fprintf('\n');



% Read data and atlas; at the same time initialize classification maps
% to atlas.
% Also construct an binary image that indicates whether a voxel should
% be classified or not based on atlas information and intensity
fprintf('Reading data and atlas      ');
data(DIM(1), DIM(2), DIM(3), nrOfChannels) = uint16(0);
atlas(DIM(1), DIM(2), DIM(3), nrOfPriors) = uint8(0);
tmpprior(DIM(1), DIM(2), DIM(3), nrOfClasses) = uint8(0);
classification(DIM(1), DIM(2), DIM(3), nrOfClasses) = uint8(0);
MRFprior(DIM(1), DIM(2), DIM(3), nrOfClasses+2) = uint8(0);
normality(DIM(1), DIM(2), DIM(3), nrOfClasses) = uint8(0);
normality(:) = 255;
lesionality(DIM(1), DIM(2), DIM(3)) = uint8(0);
garbagality(DIM(1), DIM(2), DIM(3)) = uint8(0);
playing(DIM(1), DIM(2), DIM(3)) = uint8(0);
playing(:) = 1;
for plane=1:DIM(3)
  A = spm_matrix([0 0 -plane]);
  B = A * inv(infoData(1).mat) * infoData(channel).mat;
  for channel=1:nrOfChannels
    data(:,:,plane,channel) = uint16(round( ...
	(spm_slice_vol(infoData(channel), inv(B), DIM(1:2), 1) - ...
	 minVals(channel)) / (maxVals(channel) - minVals(channel)) * 65535));
    playing(:,:,plane) = (playing(:,:,plane) & (data(:,:,plane,channel)>0));
  end
  
  atlasTmp(DIM(1), DIM(2), nrOfPriors-1) = 0;
  for prior=1:nrOfPriors-1
    B = A * inv(infoData(1).mat) * infoAtlas(prior).mat;
    atlasTmp(:,:,prior) = ...
	spm_slice_vol(infoAtlas(prior), inv(B), DIM(1:2), 0);
  end
  sumAtlasTmp = sum(atlasTmp, 3);
  overflowInd = find(sumAtlasTmp>1);
  for prior=1:nrOfPriors-1
    atlasTmp(overflowInd + (prior-1)*DIM(1)*DIM(2)) = ...
	atlasTmp(overflowInd + (prior-1)*DIM(1)*DIM(2)) ./ ...
	sumAtlasTmp(overflowInd);
  end
  for prior=1:nrOfPriors-1
    atlas(:,:,plane,prior) = uint8(255*atlasTmp(:,:,prior));
  end
  rest = 255 - sum(atlas(:,:,plane,1:nrOfPriors-1),4);
  atlas(:,:,plane,nrOfPriors) = rest;
  
  playing(:,:,plane) = (playing(:,:,plane) & (atlas(:,:,plane,nrOfPriors)<255));
  
  playingInd = find(playing(:,:,plane));
  planeInd = playingInd + (plane-1)*DIM(1)*DIM(2);
  for class=1:nrOfClasses
    classification(planeInd + (class-1)*prod(DIM)) = ...
        double(atlas(planeInd + (lkp(class)-1)*prod(DIM)))/nc(lkp(class)); 
  end
  
  fprintf('\b\b\b\b\b');
  fprintf('%-3d %%', round(100*plane/DIM(3)));
end
fprintf('\n');



% Construct a list of indices of voxels that will be sampled when
% estimating the Gaussian mixture and the bias field model.
% Sample on a subgrid that best approaches the 4x4x4 mm grid
disp('Preparing downsampling');
[dummy1 VOX dummy2 dummy3 dummy4 ORIGIN dummy6] = ...
    spm_hread(dataNames(1,:));
ideal = [4 4 4]./VOX;
downFactorTry = [floor(ideal); ceil(ideal)];
errors = 1-downFactorTry./(ones(2,1)*ideal);
[dummy ind] = min(abs(errors));
downFactor = ...
    [downFactorTry(ind(1),1) downFactorTry(ind(2),2) downFactorTry(ind(3),3)];
downFactor = max(downFactor, [1 1 1]);
if use2D
  downFactor(3) = 1;
end
downDIM = floor((DIM-1)./downFactor) + 1;

ind1 = [1:downFactor(1):DIM(1)]';    
ind2 = [1:downFactor(2):DIM(2)]';
ind3 = [1:downFactor(3):DIM(3)]';
[ind1,ind2,ind3] = ndgrid(ind1,ind2,ind3);
sampleInd = ind1 + (ind2-1)*DIM(1) + (ind3-1)*DIM(1)*DIM(2);
insideROI = find(playing(sampleInd(:)));
sampleInd = sampleInd(insideROI);

sampleData = zeros(length(sampleInd), nrOfChannels);
for channel=1:nrOfChannels
  sampleData(:,channel) = scaleFactor * ...
      log(double(data(sampleInd+(channel-1)*prod(DIM))));
end
sampleCorrectedData = zeros(length(sampleInd), nrOfChannels);
sampleClassification = zeros(length(sampleInd), nrOfClasses);
sampleNormality = zeros(length(sampleInd), nrOfClasses);
sampleLesionality = zeros(length(sampleInd),1);


% Build basis functions for the bias field model evaluated in each
% sample voxel. Calculate Q-R decomposition to orthogonalize these
% basis functions.
disp('Building basis functions for the bias field')
if ~use2D
  nrOfBasisFuncs = ...
      (maxBiasOrder+1)*(maxBiasOrder+2)/2*(maxBiasOrder+3)/3;
  sampleBasisFuncs = zeros(length(sampleInd), nrOfBasisFuncs);
  ind1 = ind1(insideROI)*2/(DIM(1)-1) - DIM(1)*2/(DIM(1)-1) + 1;
  ind2 = ind2(insideROI)*2/(DIM(2)-1) - DIM(2)*2/(DIM(2)-1) + 1;
  ind3 = ind3(insideROI)*2/(DIM(3)-1) - DIM(3)*2/(DIM(3)-1) + 1;
  basisFuncsInfo = [' '];
  ind = 1;
  for order=0:maxBiasOrder
    for xorder=0:order
      for yorder=0:order-xorder
        zorder=order-yorder-xorder;
        basisFuncsInfo = str2mat(basisFuncsInfo, ['x^' num2str(xorder) ...
              ' * y^' num2str(yorder) ' * z^' num2str(zorder)]);
        sampleBasisFuncs(:,ind) = ind1.^xorder .* ind2.^yorder .* ind3.^zorder;
        ind = ind + 1;
      end
    end
  end
  basisFuncsInfo = basisFuncsInfo(2:end,:);
  [orthoSampleBasisFuncs orthogonalizer] = qr(sampleBasisFuncs,0);
  biasCoeff = zeros(nrOfBasisFuncs, nrOfChannels);
  
  planeBasisFuncs = zeros(DIM(1)*DIM(2), maxBiasOrder+1, nrOfChannels);
  planeBiasCoeff = zeros(maxBiasOrder+1,1);
  ind1 = [1:DIM(1)]';    
  ind2 = [1:DIM(2)]';
  [ind1,ind2] = ndgrid(ind1,ind2);
  ind1 = ind1(:)*2/(DIM(1)-1) - DIM(1)*2/(DIM(1)-1) + 1;
  ind2 = ind2(:)*2/(DIM(2)-1) - DIM(2)*2/(DIM(2)-1) + 1;
  ind1powers = ones(DIM(1)*DIM(2), maxBiasOrder+1);
  ind2powers = ones(DIM(1)*DIM(2), maxBiasOrder+1);
  for order=1:maxBiasOrder
    ind1powers(:,order+1) = ind1powers(:,order) .* ind1;
    ind2powers(:,order+1) = ind2powers(:,order) .* ind2;
  end
else
  nrOfBasisFuncs = (maxBiasOrder+1)*(maxBiasOrder+2)/2;
  basisFuncs = zeros(DIM(1)*DIM(2), nrOfBasisFuncs);
  ind1 = [1:DIM(1)]';    
  ind2 = [1:DIM(2)]';
  [ind1,ind2] = ndgrid(ind1,ind2);
  ind1 = ind1(:)*2/(DIM(1)-1) - DIM(1)*2/(DIM(1)-1) + 1;
  ind2 = ind2(:)*2/(DIM(2)-1) - DIM(2)*2/(DIM(2)-1) + 1;
  basisFuncsInfo = [' '];
  ind = 1;
  for order=0:maxBiasOrder
    for xorder=0:order
      yorder=order-xorder;
      basisFuncsInfo = str2mat(basisFuncsInfo, ['x^' num2str(xorder) ...
            ' * y^' num2str(yorder)]);
      basisFuncs(:,ind) = ind1.^xorder .* ind2.^yorder;
      ind = ind + 1;
    end
  end
  basisFuncsInfo = basisFuncsInfo(2:end,:);
  biasCoeff = zeros(nrOfBasisFuncs, DIM(3), nrOfChannels);

  
  sampleBasisFuncs = zeros(length(sampleInd), nrOfBasisFuncs);
  orthoSampleBasisFuncs = zeros(length(sampleInd), nrOfBasisFuncs);
  belongsToPlane = zeros(size(sampleInd));
  for plane=1:DIM(3)
    tmp = sampleInd - (plane-1)*nrOf2DElements;
    belongsToPlaneInd = find(tmp>0 & tmp<=nrOf2DElements);
    belongsToPlane(belongsToPlaneInd) = plane;
    planeSampleInd = tmp(belongsToPlaneInd);
    sampleBasisFuncs(belongsToPlaneInd,:) = ...
        basisFuncs(planeSampleInd,:);
    
    if length(planeSampleInd)>=nrOfBasisFuncs
      [orthoSampleBasisFuncs(belongsToPlaneInd,:) orthogonalizer] = ...
          qr(sampleBasisFuncs(belongsToPlaneInd,:),0);
    end
  end
  
end


% Only use the basis functions that correspond to a bias polynomial
% order <= maxBiasOrder
biasOrder = 0;
if ~use2D
  nrOfUsedBasisFuncs = (biasOrder+1)*(biasOrder+2)/2*(biasOrder+3)/3;
else
  nrOfUsedBasisFuncs = (biasOrder+1)*(biasOrder+2)/2;
end

referenceIntensity = zeros(nrOfChannels,1);
for channel=1:nrOfChannels
  referenceIntensity(channel) = sampleData(:,channel)' * ...
      double(classification(sampleInd)) / ...
      sum(double(classification(sampleInd)));                             
end

atlasClassWeights = zeros(1, nrOfClasses);
for class=1:nrOfClasses
  atlasClassWeights(class) = sum(atlas(sampleInd + ...
      (lkp(class)-1)*nrOf3DElements));
end
atlasClassWeights = atlasClassWeights/sum(atlasClassWeights);
historyOfParams = struct(...
    'means', [], ...
    'covariances', [], ...
    'backgroundVariance', [], ...
    'classWeights', [], ...
    'biasCoeff', [], ...
    'G', [], ...
    'H', [], ...
    'relativeChangeCost', [], ...
    'atlasClassWeights', atlasClassWeights); 


[lesionConstraint, garbageConstraint, lesionHomeClass, ...
	  garbageHomeClass] = transformInputStrings(lesionConstraint, ...
						  garbageConstraint, ...
						  lesionHomeClass, ...
						  garbageHomeClass, ...
						  nrOfChannels);
relativeChangeCost = 1;
converged = 0;
forcedAtFullResolution = 0;
iteration = 0;



while (~converged | forcedAtFullResolution)


  iteration = iteration + 1;
  
  disp(['---------------']);
  if ~forcedAtFullResolution
    disp(['Iteration ' num2str(iteration)]);
    
    if (relativeChangeCost<0.003 & biasOrder<maxBiasOrder)
      disp(['Increasing bias order to ' num2str(biasOrder+1)]);
      biasOrder = biasOrder+1;
      if ~use2D
        nrOfUsedBasisFuncs = ...
            (biasOrder+1)*(biasOrder+2)/2*(biasOrder+3)/3;
      else
        nrOfUsedBasisFuncs = ...
            (biasOrder+1)*(biasOrder+2)/2;
      end
    end

    
    
    % Estimate Gaussian mixture parameters mean and covariance
    if ~use2D
      for channel=1:nrOfChannels
        sampleCorrectedData(:,channel) = sampleData(:,channel) - ...
            sampleBasisFuncs * biasCoeff(:,channel);
      end
    else
      for channel=1:nrOfChannels
        for plane=1:DIM(3)
          belongsToPlaneInd = find(belongsToPlane==plane);
          sampleCorrectedData(belongsToPlaneInd,channel) = ...
              sampleData(belongsToPlaneInd,channel) - ...
              sampleBasisFuncs(belongsToPlaneInd,:) * biasCoeff(:,plane,channel);
        end
      end
    end
    for class=1:nrOfClasses
      sampleClassification(:,class) = ...
          classification(sampleInd+(class-1)*prod(DIM));
      sampleNormality(:,class) = ...
          normality(sampleInd+(class-1)*prod(DIM));
    end
    sampleLesionality(:) = lesionality(sampleInd);

    
    % =======================================================================
    % Before re-estimating params, measure cost function with present params
    if ~(iteration==1)
      presentCost = 0;
      for class=1:nrOfClasses-1
        mahalanobis_sq = (sampleCorrectedData - ones(size(sampleInd))*means(:,class)') / ...
            sqrtm(covariances(:,:,class));
        mahalanobis_sq =  mahalanobis_sq.^2;
        if (nrOfChannels>1)
          mahalanobis_sq = sum(mahalanobis_sq, 2);
        end
        presentCost = presentCost + .5 * sampleClassification(:,class)' *  ...
              (mahalanobis_sq + log(det(covariances(:,:,class))));
      end
      if (nrOfChannels==1)
        tmp = exp(2*sampleData/scaleFactor)/backgroundVariance;
        presentCost = presentCost - sampleClassification(:,nrOfClasses)' *  ...  
            log(tmp.*exp(-tmp/2)/scaleFactor + eps);
      else
        presentCost = presentCost - sampleClassification(:,nrOfClasses)' *  ...
            log(exp(2*sum(sampleData,2)/scaleFactor)/prod(backgroundVariance) .* ...
            exp( -sum((exp(2*sampleData/scaleFactor) ./ ...
            (ones(size(sampleInd))*backgroundVariance')), 2)/2 ) / ...
            (scaleFactor^nrOfChannels) + eps);
      end
    end
    % =======================================================================
    

    disp('Estimating mixture parameters')
    for class=1:nrOfClasses
      means(:,class) = sampleCorrectedData' * ...
          sampleClassification(:,class) / ...
          sum(sampleClassification(:,class));

      %covariances(:,:,class) = (sampleCorrectedData' * ...
      %    ((sampleClassification(:,class)*ones(1,nrOfChannels)) .* ...
      %    sampleCorrectedData)) ...
      %    / sum(sampleClassification(:,class)) ...
      %    - means(:,class)*means(:,class)';

      covariances(:,:,class) = ...
          (sampleCorrectedData - ones(size(sampleInd))*means(:,class)')' * ...
          ((sampleClassification(:,class) *ones(1,nrOfChannels)) .* ...
          (sampleCorrectedData - ...
          ones(size(sampleInd))*means(:,class)')) / ...
          sum(sampleClassification(:,class));
      
      classWeights(class) = sum(sampleClassification(:,class)) / ...
          length(sampleInd)/255;
    end
    classWeights(1) = sum(sampleClassification(:,1)) / ...
          length(sampleInd)/255;
    backgroundVariance(:) = exp(2*sampleData'/scaleFactor) * ...
    sampleClassification(:,end) / sum(sampleClassification(:,end)) /2;

    
    
    for channel=1:nrOfChannels
      means(channel,:) = means(channel,:) - ...
          (means(channel,1) - referenceIntensity(channel));
    end

    if iteration==1
      % Split clusters that had the same initial classification
      for i=1:nc(end)
        covariances(:,:,3+i) = covariances(:,:,3+i)/0.8^(i-1);
        means(:,3+i) = means(:,3+i)/0.8^(i-1);
      end
      
      % Background is much darker than estimated based on all non-brain tissues
      backgroundVariance = backgroundVariance/100;
    end


    
    % Estimate bias field parameters
    if ~use2D
      disp(['Estimating bias field parameters (3D, order ' ...
            num2str(biasOrder) ')'])
    else
      disp(['Estimating bias field parameters (2D, order ' ...
            num2str(biasOrder) ')'])
    end
    invCovariances = zeros(nrOfChannels, nrOfChannels, nrOfClasses);
    for class=1:nrOfClasses
      invCovariances(:,:,class) = inv(covariances(:,:,class));
    end
    
    weights = zeros(length(sampleInd), nrOfChannels, nrOfChannels);
    for row=1:nrOfChannels
      for col=1:nrOfChannels
        weights(:,row,col) = sampleClassification(:,1:end-1) * ...
            squeeze(invCovariances(row,col,1:end-1));
      end
    end
    
    if ~use2D
      lhs = zeros(nrOfChannels*nrOfUsedBasisFuncs);
      rhs = zeros(nrOfChannels*nrOfUsedBasisFuncs, 1);
      for row=1:nrOfChannels
        tmp = squeeze(invCovariances(row,:,1:end-1));
        if nrOfChannels==1
          tmp = tmp';
        end
        rhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
            orthoSampleBasisFuncs(:,1:nrOfUsedBasisFuncs)' * ...
            sum((squeeze(weights(:,row,:)) .* sampleData - ...
            sampleClassification(:,1:end-1) * ...
            (tmp .* means(:,1:end-1))'), 2);
        for col=1:nrOfChannels
          if (col>=row)
            lhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                (col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                orthoSampleBasisFuncs(:,1:nrOfUsedBasisFuncs)' * ...
                ((weights(:,row,col) * ones(1, nrOfUsedBasisFuncs)) .* ...
                sampleBasisFuncs(:,1:nrOfUsedBasisFuncs));
          else
            lhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                (col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                lhs((col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                (row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]);    
          end
        end
      end
      solution = lhs\rhs;
      solution = reshape(solution, [nrOfUsedBasisFuncs nrOfChannels]);
      for channel=1:nrOfChannels
        biasCoeff(1:nrOfUsedBasisFuncs,channel) = solution(:,channel);
        biasCoeff(nrOfUsedBasisFuncs+1:nrOfBasisFuncs,channel) = 0;
      end
    else
      lhs = zeros(nrOfChannels*nrOfUsedBasisFuncs);
      rhs = zeros(nrOfChannels*nrOfUsedBasisFuncs,1);
      
      for plane=1:DIM(3)
        belongsToPlaneInd = find(belongsToPlane==plane);
        for row=1:nrOfChannels
          tmp = squeeze(invCovariances(row,:,1:end-1));
          if nrOfChannels==1
            tmp = tmp';
          end
          if length(belongsToPlaneInd)==1
            rhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                orthoSampleBasisFuncs(belongsToPlaneInd,1:nrOfUsedBasisFuncs)' * ...
                sum((squeeze(weights(belongsToPlaneInd,row,:))' .* ...
                sampleData(belongsToPlaneInd,:) - ...
                sampleClassification(belongsToPlaneInd,1:end-1) * ...
                (tmp .* means(:,1:end-1))'), 2);            
          else
            rhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                orthoSampleBasisFuncs(belongsToPlaneInd,1:nrOfUsedBasisFuncs)' * ...
                sum((squeeze(weights(belongsToPlaneInd,row,:)) .* ...
                sampleData(belongsToPlaneInd,:) - ...
                sampleClassification(belongsToPlaneInd,1:end-1) * ...
                (tmp .* means(:,1:end-1))'), 2);
          end
          for col=1:nrOfChannels
            if (col>=row)
              lhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                  (col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                  orthoSampleBasisFuncs(belongsToPlaneInd,1:nrOfUsedBasisFuncs)' * ...
                  ((weights(belongsToPlaneInd,row,col) * ...
                  ones(1, nrOfUsedBasisFuncs)) .* ...
                  sampleBasisFuncs(belongsToPlaneInd,1:nrOfUsedBasisFuncs));
            else
              lhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                  (col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                  lhs((col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                  (row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]);    
            end
          end
        end
        if rcond(lhs)>eps
          solution = lhs\rhs;
          solution = reshape(solution, [nrOfUsedBasisFuncs nrOfChannels]);
          for channel=1:nrOfChannels
            biasCoeff(1:nrOfUsedBasisFuncs,plane,channel) = solution(:,channel);
            biasCoeff(nrOfUsedBasisFuncs+1:nrOfBasisFuncs,plane,channel) = 0;
          end
        else
          disp(['  bias for plane ' num2str(plane) ' could not be updated']);
        end
      end
    end



    % =======================================================================
    % Measure cost function with updated params
    if ~use2D
      for channel=1:nrOfChannels
        sampleCorrectedData(:,channel) = sampleData(:,channel) - ...
            sampleBasisFuncs * biasCoeff(:,channel);
      end
    else
      for channel=1:nrOfChannels
        for plane=1:DIM(3)
          belongsToPlaneInd = find(belongsToPlane==plane);
          sampleCorrectedData(belongsToPlaneInd,channel) = ...
              sampleData(belongsToPlaneInd,channel) - ...
              sampleBasisFuncs(belongsToPlaneInd,:) * biasCoeff(:,plane,channel);
        end
      end
    end
    
    newCost = 0;
    for class=1:nrOfClasses-1
      mahalanobis_sq = (sampleCorrectedData - ones(size(sampleInd))*means(:,class)') / ...
          sqrtm(covariances(:,:,class));
      mahalanobis_sq =  mahalanobis_sq.^2;
      if (nrOfChannels>1)
        mahalanobis_sq = sum(mahalanobis_sq, 2);
      end
      newCost = newCost + .5 * sampleClassification(:,class)' *  ...
          (mahalanobis_sq + log(det(covariances(:,:,class))));
    end
    if (nrOfChannels==1)
      tmp = exp(2*sampleData/scaleFactor)/backgroundVariance;
      newCost = newCost - sampleClassification(:,nrOfClasses)' *  ...  
          log(tmp.*exp(-tmp/2)/scaleFactor + eps);
    else
      newCost = newCost - sampleClassification(:,nrOfClasses)' *  ...
          log(exp(2*sum(sampleData,2)/scaleFactor)/prod(backgroundVariance) .* ...
          exp( -sum((exp(2*sampleData/scaleFactor) ./ ...
          (ones(size(sampleInd))*backgroundVariance')), 2)/2 ) / ...
          (scaleFactor^nrOfChannels) + eps);
    end
    
    
    if (iteration==1)
      presentCost = newCost * 2;
    end
    
    % =======================================================================


    % Estimate Markov random field parameters
    if updateMRFparams
      disp('Estimating Markov random field parameters')
      [G, H] = ems_getMRFparams(classification, lkp, sampleInd, 3);
    end
    
  else
    disp('Final classification at full resolution forced')
  end
  
  
  % Classify and calculate bias field plane per plane
  fprintf('Calculating classification      ');
  if ~use2D
    ind = 1;
    planeBasisFuncs(:) = 0;
    for order=0:maxBiasOrder
      for xorder=0:order
        for yorder=0:order-xorder
          zorder=order-yorder-xorder;
          for channel=1:nrOfChannels
            planeBasisFuncs(:, zorder+1, channel) = ...
                planeBasisFuncs(:, zorder+1, channel) + ...
                biasCoeff(ind, channel) * ind1powers(:,xorder+1) ...
                .* ind2powers(:,yorder+1);
            % biasCoeff(ind, channel) * ind1.^xorder .* ind2.^yorder;
          end
          ind = ind + 1;
        end
      end
    end
  end
  
  % First store spatially varying prior in classification, and MRF
  % prior in MRFprior 
  if (useMRF | forcedAtFullResolution)
    if useMRF
      uniformAtlas(DIM(1), DIM(2), DIM(3), 4) = uint8(0);
      uniformAtlas(:) = 1;
      classification(:,:,:,lesionHomeClass) = ...
	  double(classification(:,:,:,lesionHomeClass)) + double(lesionality);
      classification(:,:,:,garbageHomeClass) = ...
	  double(classification(:,:,:,garbageHomeClass)) + double(garbagality);
      classification(:) = ems_getMRFprior(classification, playing, ...
					  uniformAtlas, lkp, G, H);
      MRFprior(:,:,:,1) = classification(:,:,:,lesionHomeClass);
      MRFprior(:,:,:,2) = classification(:,:,:,garbageHomeClass);
      MRFprior(:,:,:,3:8) = classification;
    else
      MRFprior(:) = 1;
    end
    classification(:) = ems_getMRFprior(classification, playing, ...
					atlas, lkp, zeros(nrOfPriors), ...
					zeros(nrOfPriors));
  else
    MRFprior(:) = 1;
    offset = 0;
    for class=1:nrOfClasses
      classification(sampleInd+offset) = ...
      double(atlas(sampleInd + ...
          (lkp(class)-1)*nrOf3DElements))/nc(lkp(class));
      offset = offset + nrOf3DElements;
    end
  end

  for plane=1:DIM(3)
    if (useMRF | forcedAtFullResolution)
      playingInd = find(playing(:,:,plane));
    else
      tmp = sampleInd - (plane-1)*nrOf2DElements;
      playingInd = tmp(find(tmp>0 & tmp<=nrOf2DElements)); 
    end
    planeInd = playingInd + (plane-1)*nrOf2DElements;
    
    if ~isempty(planeInd)
      if ~use2D
        ind3 = plane*2/(DIM(3)-1) - DIM(3)*2/(DIM(3)-1) + 1;
        for zorder=0:maxBiasOrder
          planeBiasCoeff(zorder+1) = ind3^zorder;
        end
      end
      
      clear planeData
      clear planeCorrectedData
      offset = 0;
      for channel=1:nrOfChannels
        planeData(:,channel) = ...
            scaleFactor * log(double(data(planeInd+offset)));
        if ~use2D
          planeCorrectedData(:,channel) = planeData(:,channel) - ...
              planeBasisFuncs(playingInd,:,channel) * planeBiasCoeff;
        else
          planeCorrectedData(:,channel) = planeData(:,channel) - ...
              basisFuncs(playingInd,:) * biasCoeff(:,plane,channel);
        end
        offset = offset + nrOf3DElements;
      end
      clear planePrior
      clear planeMRFprior
      offset = 0;
      for class=1:nrOfClasses
        planePrior(:,class) = double(classification(planeInd+offset));
        offset = offset + nrOf3DElements;
      end
      offset = 0;
      for class=1:nrOfClasses+2
        planeMRFprior(:,class) = double(MRFprior(planeInd+offset));
        offset = offset + nrOf3DElements;
      end

      clear planeClassification
      clear planeNormality

      eval(lesionConstraint)
      eval(garbageConstraint)


      clear planeAbnormalityPdf
      planeAbnormalityPdf = lambda * ...
          (planeMRFprior(:,1).*isGoodLesionIntensity + ...
          planeMRFprior(:,2).*isGoodGarbageIntensity);
      for class=1:nrOfClasses-1
        mahalanobis_sq = (planeCorrectedData - ones(size(playingInd))*means(:,class)') / ...
            sqrtm(covariances(:,:,class));
        mahalanobis_sq =  mahalanobis_sq.^2;
        if (nrOfChannels>1)
          mahalanobis_sq = sum(mahalanobis_sq, 2);
        end
        planeClassification(:,class) = planePrior(:,class) .* ...
            exp(-.5*mahalanobis_sq) / sqrt(det(covariances(:,:,class)));
        if dumping(class)
          planeNormality(:,class) = 1 - ...
              planeAbnormalityPdf ./ (planeAbnormalityPdf + ...
              planeMRFprior(:,class+2).*exp(-.5*mahalanobis_sq) + ...
              eps);
        else
          planeNormality(:,class) = 1;
        end
      end
      if (nrOfChannels==1)
        tmp = exp(2*planeData/scaleFactor)/backgroundVariance;
        planeClassification(:,nrOfClasses) = ...
            planePrior(:,nrOfClasses) .* ...
            tmp.*exp(-tmp/2)/scaleFactor;
      else
        planeClassification(:,nrOfClasses) = ...
            planePrior(:,nrOfClasses) .* ... 
            exp(2*sum(planeData,2)/scaleFactor)/prod(backgroundVariance) .* ...
            exp( -sum((exp(2*planeData/scaleFactor) ./ ...
            (ones(size(playingInd))*backgroundVariance')), 2)/2 ) / ...
            (scaleFactor^nrOfChannels);
      end
      planeNormality(:,nrOfClasses) = 1;

      
      likelihood =  sum(planeClassification,2)/255 + eps;
      lesionality(planeInd) = 255 * isGoodLesionIntensity;
      garbagality(planeInd) = 255 * isGoodGarbageIntensity;
      for class=1:nrOfClasses
        classification(planeInd + (class-1)*prod(DIM)) = ...
            planeClassification(:,class)./likelihood .* planeNormality(:,class);
        normality(planeInd + (class-1)*prod(DIM)) = ...
            255*planeNormality(:,class);
        lesionality(planeInd) =  double(lesionality(planeInd)) - ...
            double(classification(planeInd + (class-1)*prod(DIM))) .* ...
            isGoodLesionIntensity;
        garbagality(planeInd) = double(garbagality(planeInd)) - ...
            double(classification(planeInd + (class-1)*prod(DIM))) .* ...
            isGoodGarbageIntensity;
      end

      
    end
    
    fprintf('\b\b\b\b\b');
    fprintf('%-3d %%', round(100*plane/DIM(3)));
  end
  fprintf('\n');

  
  relativeChangeCost = (presentCost - ...
      newCost) / newCost;
  if ~forcedAtFullResolution
      converged = (relativeChangeCost<0.00005 & ...
          biasOrder==maxBiasOrder)  | (iteration==35);
    disp(['relativeChangeCost = ' num2str(relativeChangeCost)])
    
    historyOfParams.means(:,:,iteration) = means; 
    historyOfParams.covariances(:,:,:,iteration) = covariances;
    historyOfParams.backgroundVariance(:,:,iteration) = backgroundVariance;
    historyOfParams.classWeights(:,:,iteration) = classWeights;
    if ~use2D
      historyOfParams.biasCoeff(:,:,iteration) = biasCoeff;
    else
      historyOfParams.biasCoeff(:,:,:,iteration) = biasCoeff; 
    end
    historyOfParams.G(:,:,iteration) = G; 
    historyOfParams.H(:,:,iteration) = H;
    historyOfParams.relativeChangeCost(iteration) = relativeChangeCost;
  end
  
  if converged
    forcedAtFullResolution = ~forcedAtFullResolution;
  end


end



% Write results to file
image2world = spm_get_space(dataNames(1,:));
 
for class=1:nrOfClasses
  fileName = ...
      [spm_str_manip(dataNames(1,:),'rd') '_seg' num2str(class) '.img'];
  disp(['Writing ' fileName]);
  fid = fopen(fileName ,'w');
  spm_hwrite(fileName, DIM, VOX, 1/255, 2, 0, ORIGIN, ...
      'Segmented image');
  spm_get_space(fileName, image2world);
  fwrite(fid, classification(:,:,:,class), 'uint8');
  fclose(fid);
end

fileName = ...
    [spm_str_manip(dataNames(1,:),'rd') '_lesion.img'];
disp(['Writing ' fileName]);
fid = fopen(fileName ,'w');
spm_hwrite(fileName, DIM, VOX, 1/255, 2, 0, ORIGIN, ...
    'Segmented image');
spm_get_space(fileName, image2world);
fwrite(fid, lesionality, 'uint8');
fclose(fid);

fileName = ...
    [spm_str_manip(dataNames(1,:),'rd') '_garbage.img'];
disp(['Writing ' fileName]);
fid = fopen(fileName ,'w');
spm_hwrite(fileName, DIM, VOX, 1/255, 2, 0, ORIGIN, ...
    'Segmented image');
spm_get_space(fileName, image2world);
fwrite(fid, garbagality, 'uint8');
fclose(fid);


planeBias = zeros(DIM(1), DIM(2));
for channel=1:nrOfChannels
  fileName = [spm_str_manip(dataNames(channel,:),'rd') '_bias.img'];
  disp(['Writing ' fileName]);
  fid = fopen(fileName, 'w');
  spm_hwrite(fileName, DIM, VOX, 10/65535, 4, 0, ORIGIN, ...
      'Estimated bias');
  spm_get_space(fileName, image2world);

  for plane=1:DIM(3)
    planeBias(:) = 0;
    playingInd = find(playing(:,:,plane));
    if ~use2D
      ind3 = plane*2/(DIM(3)-1) - DIM(3)*2/(DIM(3)-1) + 1;
      for zorder=0:maxBiasOrder
        planeBiasCoeff(zorder+1) = ind3^zorder;
      end
      planeBias(playingInd) = ...
          exp(planeBasisFuncs(playingInd,:,channel) * ...
          planeBiasCoeff / scaleFactor);
    else
      planeBias(playingInd) = ...
          exp(basisFuncs(playingInd,:) * ...
          biasCoeff(:,plane,channel) / scaleFactor);     
    end
    planeBias = round(planeBias / 10 * 65535);
    tmp = find(planeBias > 65535);
    planeBias(tmp) = 65535;
    fwrite(fid, planeBias, 'uint16');
  end
  fclose(fid);
end  


matName = [spm_str_manip(dataNames(1,:),'rd') '_inputArguments.mat'];
disp(['Writing ' matName])
eval(['save ' matName ' inputArguments'])

matName = [spm_str_manip(dataNames(1,:),'rd') '_historyOfParams.mat'];
disp(['Writing ' matName])
eval(['save ' matName ' historyOfParams'])

disp('done');
return    


  









% ===========================================================================
% ===========================================================================
function uiForems_lesions()

SPMid = spm('FnBanner',mfilename,'2.9');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Segment');

n     = spm_input('number of subjects',1,'e',1);
if n < 1,
  spm_figure('Clear','Interactive');
  return;
end;


for i = 1:n,
  dataNames = spm_get(Inf,'.img',...
      ['Select MRI(s) for subject ' num2str(i)]);
  eval(['dataNames' num2str(i) ' = dataNames;']);
end;

maxBiasOrder = spm_input('Bias polynomial order?',1,'w',4);

use2D = spm_input('What kind of polynomial?', 2, 'b', ['2D'; '3D'], ...
    [1; 0]);

mahalanobis = spm_input('Mahalanobis threshold?',3,'r',3);

useMRF = spm_input('Use Markov random field?', 4, 'y/n', [1; 0]);

if useMRF
  updateMRFparams = 1;
else
  updateMRFparams = 0;
end

G = zeros(4,4);
H = zeros(4,4);

lesionConstraint = spm_input('Intensity constraint on lesions', 5, ...
			     's', ' ');


if ~useMRF
  lesionHomeClass = 'wm';  % Don't care
  pos = 6;
else
  lesionHomeClass = spm_input('To which tissue do lesions belong?', ...
			      6, 'b', str2mat('wm', 'gm', 'csf'), ...
			      str2mat('wm', 'gm', 'csf'));
  pos = 7;
end


tmp = spm_input('Look for other outliers as well?', pos, 'y/n', [1; 0]);
if tmp
  garbageConstraint = spm_input('Intensity constraint on other outliers', ...
				pos+1, 's', ' ');
  if ~useMRF
    garbageHomeClass = 'csf'; % Don't care
  else
    garbageHomeClass = spm_input(...
	'To which tissue do other outliers belong?', ...
	pos+2, 'b', str2mat('wm', 'gm', 'csf'), ...
	str2mat('wm', 'gm', 'csf'));
  end
else
  garbageConstraint = 'i1>Inf';
  garbageHomeClass = 'csf' % Don't care
end


  
spm('Pointer','Watch');
for i = 1:n,
  spm('FigName',['Segment: working on subj ' num2str(i)],Finter,CmdLine);
  fprintf('\Segmenting Subject %d: ', i);
  
  eval(['dataNames = dataNames' num2str(i) ';']);
  if (size(dataNames,1)~=0) 
    ems_lesions(dataNames, use2D, maxBiasOrder, useMRF, ...
        updateMRFparams, G, H, mahalanobis, lesionConstraint, ...
        lesionHomeClass, garbageConstraint, garbageHomeClass);
  end
end;

fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
spm_figure('Clear',spm_figure('FindWin','Interactive'));
spm('FigName','Segment: done',Finter,CmdLine);
spm('Pointer');
return;









% ===========================================================================
% ===========================================================================
function checkConsistency(dataNames, use2D, maxBiasOrder, useMRF, ...
			  updateMRFparams, G, H, mahalanobis, ...
			  lesionConstraint, lesionHomeClass, ...
			  garbageConstraint, garbageHomeClass)

if ~isstr(dataNames)
  error('dataNames must be string input')
end
nrOfChannels = size(dataNames,1);
for channel=1:nrOfChannels
  if (exist(deblank(dataNames(channel,:)), 'file')~=2)
    error([deblank(dataNames(channel,:)) ' not found'])
  end
end

if 0
  [refDIM refVOX refSCALE refTYPE refOFFSET refORIGIN] = ...
      spm_hread(deblank(dataNames(1,:)));
  refImage2world = spm_get_space(dataNames(1,:));
  
  for channel=2:nrOfChannels
    [DIM VOX SCALE TYPE OFFSET ORIGIN] = ...
	spm_hread(deblank(dataNames(channel,:)));
    if (~all(DIM==refDIM) | ~all(VOX==refVOX))
      error('All input files must have the same spatial grid')
    end
    image2world = spm_get_space(deblank(dataNames(channel,:)));
    if ~all(image2world(:)==refImage2world(:))
      error('All input files must have same image-to-world transformation');
    end
  end
end
  
if (~all(size(use2D)==[1 1]) | isstr(use2D))
  error('use2D must be scalar');
end
if (~all(size(maxBiasOrder)==[1 1]) | isstr(maxBiasOrder))
  error('maxBiasOrder must be scalar');
end
if (~all(size(useMRF)==[1 1]) | isstr(useMRF))
  error('useMRF must be scalar');
end
if (~all(size(updateMRFparams)==[1 1]) | isstr(updateMRFparams))
  error('updateMRFparams must be scalar');
end
if (~all(size(G)==[4 4]) | isstr(G) | ~isreal(G))
  error('G must be a [4x4] real matrix');
end
if (~all(size(H)==[4 4]) | isstr(H) | ~isreal(H))
  error('H must be a [4x4] real matrix');
end
if (~all(size(mahalanobis)==[1 1]) | isstr(mahalanobis))
  error('mahalanobis must be scalar');
end


nrOfChannels = size(dataNames,1);
[lesionConstraint, garbageConstraint, lesionHomeClass, ...
	  garbageHomeClass] = transformInputStrings(lesionConstraint, ...
						  garbageConstraint, ...
						  lesionHomeClass, ...
						  garbageHomeClass, ...
						  nrOfChannels);
return



% ===========================================================================
% ===========================================================================
function prettyEcho(dataNames, use2D, maxBiasOrder, useMRF, ...
		    updateMRFparams, G, H, mahalanobis, lesionConstraint, ...
		    lesionHomeClass, garbageConstraint, garbageHomeClass)

fprintf('\n\n');
disp('oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
disp('                                        _                     ');
disp('                                       |_|                    ');
disp('    ___  __  __  ___     _    ___  ___  _  ___  __  _  ___    ');
disp('   | __)|  \/  |/ __)   | |  | __)/ __)| |(   )|  \| |/ __)   ');
disp('   | __||      |\__ \   | |_ | __|\__ \| || | ||     |\__ \   ');
disp('   |___)|_|\/|_|(___/___|___)|___)(___/|_|(___)|_|\__|(___/   ');
fprintf('\n\n');
disp('               Koen Van Leemput       August 17, 2001         ');
fprintf('\n');
disp('                      koen.vanleemput@hut.fi                  ');
fprintf('\n');
disp('oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
fprintf('\n\n');


fprintf('      dataNames = [\n');
for i=1:size(dataNames,1)
  fprintf('        ''%s'';\n', dataNames(i,:));
end
fprintf('                  ]\n');

fprintf('      use2D = %d\n', use2D);
fprintf('      maxBiasOrder = %d\n', maxBiasOrder);
fprintf('      useMRF = %d\n', useMRF);
if useMRF
  fprintf('      updateMRFparams = %d\n', updateMRFparams);
  if ~updateMRFparams
    fprintf('      G = [\n');
    for i=1:size(G,1)
      fprintf('        %f ', G(i,:));
      fprintf('      \n');
    end
    fprintf('          ]\n');
    
    fprintf('      H = [\n');
    for i=1:size(H,1)
      fprintf('        %f ', H(i,:));
      fprintf('      \n');
    end
    fprintf('          ]\n');
  end
end
fprintf('      mahalanobis = %f\n', mahalanobis);
fprintf('      lesionConstraint = ''%s''\n', lesionConstraint);
fprintf('      lesionHomeClass = ''%s''\n', lesionHomeClass);
fprintf('      garbageConstraint = ''%s''\n', garbageConstraint);
fprintf('      garbageHomeClass = ''%s''\n', garbageHomeClass);
fprintf('\n\n');
return





% ===========================================================================
% ===========================================================================
function [lesionConstraint, garbageConstraint, lesionHomeClass, ...
	  garbageHomeClass] = transformInputStrings(lesionConstraint, ...
						  garbageConstraint, ...
						  lesionHomeClass, ...
						  garbageHomeClass, ...
						  nrOfChannels)
%
% Works only for nr of channels < 10
%  


% First make strings that implement the intensity constraints
if isempty(deblank(lesionConstraint))
  lesionConstraint = 'i1>-Inf';
end
lesionConstraintTest = replaceCodes(['isGoodLesionIntensity = '  ...
		    lesionConstraint ';']);

if isempty(deblank(garbageConstraint))
  garbageConstraint = 'i1>-Inf';
end
garbageConstraintTest = replaceCodes(['isGoodGarbageIntensity = '  ...
		    garbageConstraint ';']);


% Check if the intensity constraints can be evaluated
means = zeros(nrOfChannels, 3);
planeCorrectedData = zeros(2, nrOfChannels);
try 
  eval(lesionConstraintTest)
catch
  error(['Invalid lesionConstraint: "' lesionConstraint '"'])
end
lesionConstraint = lesionConstraintTest;

try 
  eval(garbageConstraintTest)
catch
  error(['Invalid garbageConstraint "' garbageConstraint '"'])
end
garbageConstraint = garbageConstraintTest;


 
% Now convert strings describing home classes into numbers
lesionHomeClassTest = strmatch(deblank(lesionHomeClass), ...
			       strvcat('wm', 'gm', 'csf'), 'exact'); 
if isempty(lesionHomeClassTest) 
  error(['Unknown lesionHomeClass: "' lesionHomeClass '"'])
else
  lesionHomeClass = lesionHomeClassTest;
end
garbageHomeClassTest = strmatch(deblank(garbageHomeClass), ...
				strvcat('wm', 'gm', 'csf'), 'exact'); 
if isempty(garbageHomeClassTest) 
  error(['Unknown garbageHomeClass: "' garbageHomeClass '"'])
else
  garbageHomeClass = garbageHomeClassTest;
end

return





% ===========================================================================
% ===========================================================================
function constraintString = replaceCodes(constraintString)
%

% Replace i1, i2 etc by planeCorrectedData(:,1), planeCorrectedData(:,2) etc
for channel=1:9
  ind = findstr(constraintString, ['i' num2str(channel)]);
  if ~isempty(ind)
    ind = ind(1);
    constraintString = [constraintString(1:ind-1) ...
			'planeCorrectedData(:,' num2str(channel) ')' ...
                        constraintString(ind+2:end)];
  end
end


% Replace wm1, wm2 etc by means(1,1), means(2,1) etc
for channel=1:9
  inds = findstr(constraintString, ['wm' num2str(channel)]);
  for i=1:length(inds)
    ind = inds(i);
    constraintString = [constraintString(1:ind-1) ...
			'means(' num2str(channel) ', 1)' ...
                        constraintString(ind+3:end)];
    inds = inds + 8;
  end
end

% Replace gm1, gm2 etc by means(1,2), means(2,2) etc
for channel=1:9
  inds = findstr(constraintString, ['gm' num2str(channel)]);
  for i=1:length(inds)
    ind = inds(i);
    constraintString = [constraintString(1:ind-1) ...
			'means(' num2str(channel) ', 2)' ...
                        constraintString(ind+3:end)];
    inds = inds + 8;
  end
end

% Replace csf1, csf2 etc by means(1,3), means(2,3) etc
for channel=1:9
  inds = findstr(constraintString, ['csf' num2str(channel)]);
  for i=1:length(inds)
    ind = inds(i);
    constraintString = [constraintString(1:ind-1) ...
			'means(' num2str(channel) ', 3)' ...
                        constraintString(ind+4:end)];
    inds = inds + 7;
  end
end
return





