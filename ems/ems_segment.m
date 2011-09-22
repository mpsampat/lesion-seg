function [inputArguments, historyOfParams] = ems_segment(dataNames, use2D, maxBiasOrder, useMRF, updateMRFparams, G, H)
%
% FORMAT [inputArguments, historyOfParams] = ems_segment(dataNames, ...
%               use2D, maxBiasOrder, useMRF, updateMRFparams, G, H)
%
% Model-based tissue classification of normal brain MR images, using
% an EM-algorithm that iteratively interleaves classification with
% estimation of tissue-specific intensity distribution parameters,
% polynomial bias field estimation, and estimation of Markov random
% field parameters, as described in:
% 
%   K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens, Automated
%   model-based tissue classification of MR images of the brain , IEEE
%   transactions on medical imaging, vol. 18, no. 10, pp. 897-908,
%   October 1999
%
% and
% 
%   K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens, Automated
%   model-based bias field correction of MR images of the brain , IEEE
%   transactions on medical imaging, vol. 18, no. 10, pp. 885-896,
%   October 1999
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
%
% - inputArguments: contains the input arguments with which the
%   program was called
% - historyOfParams: contains the estimated model parameters at every
%   iteration.
%
% General remarks:
% - if no bias correction is desired, select use2D=0 and maxBiasOrder=0.
% 
% ------------------------------------------------------------------------
% ems_segment.m    Koen Van Leemput - August 17, 2001



% Get input data, check and correct if necessary
if (nargin==6)
  H = zeros(4,4);
elseif (nargin==5)
  G = zeros(4,4);
  H = zeros(4,4); 
elseif (nargin==4)
  if useMRF
    updateMRFparams = 1;
  else
    updateMRFparams = 0;
  end
  G = zeros(4,4);
  H = zeros(4,4);
elseif (nargin==3)
  useMRF = 0;
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
elseif (nargin==2)
  useMRF = 0;
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
  maxBiasOrder = 4;
elseif (nargin==1)
  useMRF = 0;
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
  maxBiasOrder = 4;
  use2D = 0;
elseif (nargin==0)
  uiForEms_segment
  return
end

if ~useMRF
  updateMRFparams = 0;
  G = zeros(4,4);
  H = zeros(4,4);
end
checkConsistency(dataNames, use2D, maxBiasOrder, useMRF, updateMRFparams, G, H)
prettyEcho(dataNames, use2D, maxBiasOrder, useMRF, updateMRFparams, G, H)



% Some initialization and allocation of variables
scaleFactor=100;

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
    'H', H);



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
classification(DIM(1), DIM(2), DIM(3), nrOfClasses) = uint8(0);
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
      (lkp(class)-1)*nrOf3DElements))/nc(lkp(class));
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
    'logLikelihood', [], ...
    'atlasClassWeights', atlasClassWeights); 
    
logLikelihood = 1/eps;
relativeChangeLogLikelihood = 1;
converged = 0;
forcedAtFullResolution = 0;
iteration = 0;


while (~converged | forcedAtFullResolution)

  iteration = iteration + 1;
  previousLogLikelihood = logLikelihood;
  
  disp(['---------------']);
  if ~forcedAtFullResolution
    disp(['Iteration ' num2str(iteration)]);
    
    if (relativeChangeLogLikelihood<0.003 & biasOrder<maxBiasOrder)
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
    disp('Estimating mixture parameters')
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
    end
    
    for class=1:nrOfClasses
      means(:,class) = sampleCorrectedData' * ...
          sampleClassification(:,class) / ...
          sum(sampleClassification(:,class));
      covariances(:,:,class) = (sampleCorrectedData' * ...
          ((sampleClassification(:,class)*ones(1,nrOfChannels)) .* ...
          sampleCorrectedData)) ...
          / sum(sampleClassification(:,class)) ...
          - means(:,class)*means(:,class)';
      classWeights(class) = sum(sampleClassification(:,class)) / ...
          length(sampleInd)/255;
      
    end
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
          rhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
              orthoSampleBasisFuncs(belongsToPlaneInd,1:nrOfUsedBasisFuncs)' * ...
              sum((squeeze(weights(belongsToPlaneInd,row,:)) .* ...
              sampleData(belongsToPlaneInd,:) - ...
              sampleClassification(belongsToPlaneInd,1:end-1) * ...
              (tmp .* means(:,1:end-1))'), 2);
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
    
    
    
    % Estimate Markov random field parameters
    if updateMRFparams
      disp('Estimating Markov random field parameters')
      [G, H] = ems_getMRFparams(classification, lkp, sampleInd, 1);
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
  
  % First store spatially varying prior in classification 
  if (useMRF | forcedAtFullResolution)
    if useMRF
      classification(:) = ...
          ems_getMRFprior(classification, playing, atlas, lkp, G, H);
    else
      classification(:) = ...
          ems_getMRFprior(classification, playing, atlas, lkp, 0*G, 0*H);
    end
  else
    offset = 0;
    for class=1:nrOfClasses
      classification(sampleInd+offset) = ...
      double(atlas(sampleInd + ...
          (lkp(class)-1)*nrOf3DElements))/nc(lkp(class));
      offset = offset + nrOf3DElements;
    end
  end
  
  logLikelihood = 0;
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
      offset = 0;
      for class=1:nrOfClasses
        planePrior(:,class) = double(classification(planeInd+offset));
        offset = offset + nrOf3DElements;
      end
      
      clear planeClassification
      for class=1:nrOfClasses-1
        mahalanobis_sq = (planeCorrectedData - ones(size(playingInd))*means(:,class)') / ...
            sqrtm(covariances(:,:,class));
        mahalanobis_sq =  mahalanobis_sq.^2;
        if (nrOfChannels>1)
          mahalanobis_sq = sum(mahalanobis_sq, 2);
        end
        planeClassification(:,class) = ...
            planePrior(:,class) .* ...
            exp(-.5*mahalanobis_sq) / sqrt(det(covariances(:,:,class)));
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
      
      likelihood =  sum(planeClassification,2)/255 + eps;
      for class=1:nrOfClasses
        classification(planeInd + (class-1)*prod(DIM)) = ...
            planeClassification(:,class)./likelihood;
      end
      
      tmp = sampleInd - (plane-1)*nrOf2DElements;
      planeSampleInd = tmp(find(tmp>0 & tmp<=nrOf2DElements)); 
      tmp = zeros(DIM(1),DIM(2));
      tmp(playingInd) = likelihood;
      logLikelihood = logLikelihood - sum(log(tmp(planeSampleInd)));
    end
    
    fprintf('\b\b\b\b\b');
    fprintf('%-3d %%', round(100*plane/DIM(3)));
  end
  fprintf('\n');
  
  % Convergence detected?
  relativeChangeLogLikelihood = (previousLogLikelihood - ...
      logLikelihood) / logLikelihood;
  if ~forcedAtFullResolution
    converged = (relativeChangeLogLikelihood<0.0001 & ...
        biasOrder==maxBiasOrder) | (iteration==35);
    disp(['-logLikelihood = ' num2str(logLikelihood)])

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
    historyOfParams.logLikelihood(iteration) = logLikelihood;

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
function uiForEms_segment()

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

maxBiasOrder = spm_input('Bias polynomial order?',2,'w',4);

use2D = spm_input('What kind of polynomial?', 3, 'b', ['2D'; '3D'], ...
    [1; 0]);

useMRF = spm_input('Use Markov random field?', 4, 'y/n', [1; 0]);

if useMRF
  updateMRFparams = 1;
else
  updateMRFparams = 0;
end

G = zeros(4,4);
H = zeros(4,4);

  
spm('Pointer','Watch');
for i = 1:n,
  spm('FigName',['Segment: working on subj ' num2str(i)],Finter,CmdLine);
  fprintf('\Segmenting Subject %d: ', i);
  
  eval(['dataNames = dataNames' num2str(i) ';']);
  if (size(dataNames,1)~=0) 
    ems_segment(dataNames, use2D, maxBiasOrder, useMRF, ...
        updateMRFparams, G, H);
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
    updateMRFparams, G, H)

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


return



% ===========================================================================
% ===========================================================================
function prettyEcho(dataNames, use2D, maxBiasOrder, useMRF, updateMRFparams, G, H)

fprintf('\n\n');
disp('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
fprintf('\n\n');
disp('    ___  __  __  ___     ___  ___  ____  __  __  ___  __  _  _____   ')
disp('   | __)|  \/  |/ __)   / __)| __)/  __)|  \/  || __)|  \| |(_   _)  ')
disp('   | __||      |\__ \   \__ \| __|| (_ ||      || __||     |  | |    ')
disp('   |___)|_|\/|_|(___/___(___/|___)\____||_|\/|_||___)|_|\__|  |_|    ') 
fprintf('\n\n');
disp('               Koen Van Leemput       August 17, 2001                ');
fprintf('\n');
disp('                      koen.vanleemput@hut.fi                         ');
fprintf('\n');
disp('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
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
fprintf('\n\n');
return
