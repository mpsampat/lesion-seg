function ems_mireg(objectName, targetName, others, free, startingParams)
% 
% FORMAT ems_mireg(objectName, targetName, others, free, ...
% 		 startingParams)
% 
% Affine registration of multi-modal images based on maximization
% of mutual information, as described in: 
%
%   F. Maes, A. Collignon, D. Vandermeulen, G. Marchal, P. Suetens,
%   Multimodality image registration by maximization of mutual
%   information , IEEE transactions on Medical Imaging, vol. 16, no. 2,
%   pp. 187-198, April 1997 (in December 2000 recognized by the IEEE as
%   the journal's most frequently cited paper published in 1997)
%
% The code is partly copied from spm_mireg.m (John Ashburner)
%
% - objectName: object image (the image to reposition).
% - targetName: target image (the image to act as template).
% - others    : other images that will undergo the same transformation
%               as objectName  
% - free      : one for free parameters, zero otherwise [1x12]
% - startingParams: starting parameters
%
% ------------------------------------------------------------------------
% ems_mireg.m    Koen Van Leemput - August 23, 2001


% Hard-coded optimization parameters
global OPTIMPARAMS
OPTIMPARAMS = struct('maxBracketStep', 1.00e+01, ...
		     'absoluteToleranceBrent', 1.00e-03, ...
		     'fractionalToleranceBrent', 1.0e-2, ...
		     'maxIterationsBrent', 30, ...
		     'fractionalTolerancePowell', 1e-5, ...
		     'maxIterationsPowell', 10, ...
		     'acceptNewDirection', 0);


% Check input arguments
if (nargin==4)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
elseif (nargin==3)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  free = ones(1,12);
elseif (nargin==2)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  free = ones(1,12);
  others = '';
elseif (nargin==1)
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
  free = ones(1,12);
  others = '';

  global SWD
  targetName = fullfile(SWD, 'templates/T1.img');
elseif (nargin==0)
  uiForEms_mireg
  return
end

if ((size(objectName,1)~=1) | (size(targetName,1)~=1))
  error('objectName and targetName must be single images');
end

prettyEcho(objectName, targetName, others, free, startingParams)


% Load the data
targetInfo = spm_vol(targetName);
objectInfo = spm_vol(objectName);
if ~isfield(targetInfo, 'uint8')
  disp(['Loading and converting ' deblank(targetName) ' into uint8'])
  targetInfo.uint8 = loaduint8(targetInfo);
end
if ~isfield(objectInfo, 'uint8')
  disp(['Loading and converting ' deblank(objectName) ' into uint8'])
  objectInfo.uint8 = loaduint8(objectInfo);
end


% Rotation around z-axis is best performed before z-translation
order = zeros(1,12);
free = [free(1:2) free(6) free(3:5) free(7:12)]; 
order = cumsum(free).*free;
order = [order(1:2) order(4:6) order(3) order(7:12)];

% Set initial search directions
x  = startingParams(:);
origXi = eye(12, 12);
origSteps = 5*[1 1 1 0.0175 0.0175 0.0175 ...
	   .02 .02 .02 .02 .02 .02];
xi = zeros(12, length(find(free)));
steps = zeros(1, length(find(free)));
for paramNr=1:12
  if order(paramNr)
    xi(:,order(paramNr)) = origXi(:,paramNr);
    steps(:, order(paramNr)) = origSteps(:,paramNr);
  end
end




% A 2-resolution strategy is applied, where the downsampling best
% approximates a 4x4x4 mm^3 grid  
[dummy1 VOX dummy2 dummy3 dummy4 dummy5 dummy6] = spm_hread(targetName);
ideal = [4 4 4]./VOX;
downFactorTry = [floor(ideal); ceil(ideal)];
errors = 1-downFactorTry./(ones(2,1)*ideal);
[dummy ind] = min(abs(errors));
downFactor = ...
    [downFactorTry(ind(1),1) downFactorTry(ind(2),2) downFactorTry(ind(3),3)];
downFactor = max(downFactor, [1 1 1]);


% Search the transformation that maximizes MI 
tic
fprintf('Initial position:\n')
printPosition(x)
[x, fval, dummy] = powell(x, xi, steps, ...
			  targetInfo, objectInfo, downFactor);
[x, fval, dummy] = powell(x, xi, steps/10, ...
			  targetInfo, objectInfo, [1 1 1]);
fprintf('Final position:\n')
printPosition(x)
t = toc;
disp(['Elapsed time: ' num2str(t) ' seconds'])


% Show some graphical output
display_results(targetInfo, objectInfo, x);


% Finally, apply the estimated transformation
M = inv(spm_matrix(x));
if isempty(others)
  Images = objectName;
else
  Images   = str2mat(objectName,others);
end
for i=1:size(Images,1)
  MM = spm_get_space(deblank(Images(i,:)));
  spm_get_space(deblank(Images(i,:)), M*MM);
end

return








% ===========================================================================
% ===========================================================================
function [p, fret, xi] = powell(p, xi, steps, targetInfo, objectInfo, s)
% Powell optimisation method - taken from Numerical Recipes (p. 417) and
% modified slightly.

global nrOfEvaluations
nrOfEvaluations=0;
p=p(:);

global OPTIMPARAMS
ITMAX = OPTIMPARAMS.maxIterationsPowell;
ftol = OPTIMPARAMS.fractionalTolerancePowell;
acceptNewDirection = OPTIMPARAMS.acceptNewDirection;

fret  = getMI(p, targetInfo, objectInfo, s);
pt    = p;
for iter=1:ITMAX,
	fp   = fret;
	ibig = 0;
	del  = 0.0;
	for i=1:size(xi, 2),
		fptt = fret;
		fprintf(['Downsampling [%d %d %d], iteration %d, ' ...
			 'search direction:\n'], s(1), s(2), s(3), iter)
		[p,xit,fret] = linmin(p, xi(:,i), steps(i), targetInfo, objectInfo, s);
		if abs(fptt-fret) > del,
			del  = abs(fptt-fret);
			ibig = i;
		end;
	end;
	steps = steps/2;
	if 2.0*abs(fp-fret) <= ftol*(abs(fp)+abs(fret)),
		return;
	end;
	ptt  = 2.0*p-pt;
	xit  = p-pt;
	pt   = p;
	fptt = getMI(ptt, targetInfo, objectInfo, s);
	if fptt < fp,
		t = 2.0*(fp-2.0*fret+fptt)*(fp-fret-del).^2-del*(fp-fptt).^2;
		if t < 0.0,
		        fprintf(['Downsampling [%d %d %d], iteration %d, ' ...
			   'search direction:\n'], s(1), s(2), s(3), iter)
			[p,xit,fret] = ...
			    linmin(p, xit, 1, targetInfo, objectInfo, s);
			if acceptNewDirection
			  xi(:,ibig)   = xi(:,end);
			  steps(:,ibig) = steps(:,end);
			  xi(:,end)    = xit;
			  steps(:,end) = 1;
			end
		end;
	end;
end;
%warning('Too many iterations in routine POWELL');
return;




% ===========================================================================
% ===========================================================================
function [p,xi,fret] = linmin(p,xi,step, targetInfo, objectInfo, s)
% Code based on Numerical Recipes in C (p. 419)

printPosition(xi)

global OPTIMPARAMS
global lnm
lnm = struct('pcom', p, 'xicom', xi, ...
	     'targetInfo', targetInfo, 'objectInfo', objectInfo, ...
	     's', s);
ax    = 0.0;
xx    = step;
[ax,xx,bx,fa,fx,fb] = mnbrak(ax,xx);
[fret,xmin] = brent(ax, xx, bx, fx);
xi    = xi * xmin;
p     = p + xi;
global nrOfEvaluations
fprintf(' Updated position:        %15.10f     (%d evaluations)\n', ...
	-fret, nrOfEvaluations)
printPosition(p)

return;





% ===========================================================================
% ===========================================================================
function [ax,bx,cx,fa,fb,fc] = mnbrak(ax,bx)
% Code based on Numerical Recipes in C (p. 400)

global OPTIMPARAMS
global lnm
pcom = lnm.pcom;
xicom = lnm.xicom;
targetInfo = lnm.targetInfo;
objectInfo = lnm.objectInfo;
s = lnm.s;


GOLD   = 1.618034;
GLIMIT = OPTIMPARAMS.maxBracketStep;
TINY   = 1.0e-20;

fprintf(' Bracket\n')

fa = getMI(pcom + ax.*xicom, targetInfo, objectInfo, s);
fb = getMI(pcom + bx.*xicom, targetInfo, objectInfo, s);

if fb > fa
	dum = ax; ax = bx; bx = dum;
	dum = fb; fb = fa; fa = dum;
end;
cx = bx+GOLD*(bx-ax);
fc = getMI(pcom + cx.*xicom, targetInfo, objectInfo, s);
printInterval(ax, cx, bx, -fb)
while fb > fc,
	r    = (bx-ax)*(fb-fc);
	q    = (bx-cx)*(fb-fa);
	u    = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*max(abs(q-r),TINY)*sign(q-r));
	ulim = bx+GLIMIT*(cx-bx);
	if (bx-u)*(u-cx) > 0.0,
		fu = getMI(pcom + u.*xicom, targetInfo, objectInfo, s);
		
		if fu < fc,
			ax = bx; bx =  u;
			fa = fb; fb = fu;
			printInterval(ax, cx, bx, -fb)
			return;
		elseif fu > fb,
			cx = u;
			fc = fu;
			printInterval(ax, cx, bx, -fb)
			return;
		end;
		u  = cx+GOLD*(cx-bx);
		fu = getMI(pcom + u.*xicom, targetInfo, objectInfo, s);	
	elseif (cx-u)*(u-ulim) > 0.0
		fu = getMI(pcom + u.*xicom, targetInfo, ...
					 objectInfo, s);
		if fu < fc,
			bx = cx; cx = u; u = cx+GOLD*(cx-bx);
			fb = fc; fc = fu; fu = getMI(pcom + u.*xicom, targetInfo, objectInfo, s);	
			printInterval(ax, cx, bx, -fb)
		end;
	elseif (u-ulim)*(ulim-cx) >= 0.0,
		u  = ulim;
		fu = getMI(pcom + u.*xicom, targetInfo, objectInfo, s);	
	else,
		u  = cx+GOLD*(cx-bx);
		fu = getMI(pcom + u.*xicom, targetInfo, objectInfo, s);	
	end;
	ax = bx; bx = cx; cx = u;
	fa = fb; fb = fc; fc = fu;
	printInterval(ax, cx, bx, -fb)
end;
return;





% ===========================================================================
% ===========================================================================
function [fx, x] = brent(ax,bx,cx,fx)
% Code based on Numerical Recipes in C (p. 404)

global OPTIMPARAMS
ITMAX = OPTIMPARAMS.maxIterationsBrent;
tol3 = OPTIMPARAMS.absoluteToleranceBrent;
tol = OPTIMPARAMS.fractionalToleranceBrent;

global lnm
pcom = lnm.pcom;
xicom = lnm.xicom;
targetInfo = lnm.targetInfo;
objectInfo = lnm.objectInfo;
s = lnm.s;

CGOLD = 0.3819660; % 1-(1-sqrt(5))/2
e = 0.0;
a = min(ax,cx);
b = max(ax,cx);
x = bx; w = bx; v = bx;
fw = fx;
fv = fx;
fprintf(' Brent\n')
printInterval(a, b, x, -fx)
for iter=1:ITMAX,
	xm   = 0.5*(a+b);
	tol1 = tol*abs(x)+eps;
	tol2 = 2.0*tol1;
	if ((abs(x-xm) + .5*(b-a) <= tol2) | ...
	    (abs(x-xm) + .5*(b-a) <= tol3))
		return;
	end;
	if abs(e) > tol1,
		r     = (x-w)*(fx-fv);
		q     = (x-v)*(fx-fw);
		p     = (x-v)*q-(x-w)*r;
		q     = 2.0*(q-r);
		if q > 0.0, p = -p; end;
		q     = abs(q);
		etemp = e;
		e     = d;
		if abs(p) >= abs(0.5*q*etemp) | p <= q*(a-x) | p >= q*(b-x),
			if x >= xm, e = a-x; else, e = b-x; end;
			d = CGOLD*(e);
		else,
			d = p/q;
			u = x+d;
			if u-a < tol2 | b-u < tol2,
				d = tol1*sign(xm-x);
			end;
		end;
	else,
		if x>=xm, e = a-x; else, e = b-x; end;
		d = CGOLD*e;
	end;
	if abs(d) >= tol1, u = x+d; else, u = x+tol1*sign(d); end;
	fu = getMI(pcom + u.*xicom, targetInfo, objectInfo, s);
	if fu <= fx,
		if u >= x, a=x; else, b=x; end;
		v =  w;  w =  x;  x =  u;
		fv = fw; fw = fx; fx = fu;
		printInterval(a, b, x, -fx)
	else,
		if u < x, a=u; else, b=u; end;
		printInterval(a, b, x, -fx)
		if fu <= fw | w == x,
			v  = w;  w =  u;
			fv = fw; fw = fu;
		elseif fu <= fv | v == x | v == w,
			v =  u;
			fv = fu;
		end;
	end;
end;
%warning('Too many iterations in BRENT');
return;




% ===========================================================================
% ===========================================================================
function o = getMI(x, targetInfo, objectInfo, s)
% The function that is minimised.

% Calculate 2-D histogram 
H = spm_hist2(targetInfo.uint8, objectInfo.uint8, ...
	      objectInfo.mat\spm_matrix(x(:)')*targetInfo.mat, s);

% Compute mutual information, discarding zero bins to reduce weight of
% background voxels on parameter estimation
H = H(2:end, 2:end);
H  = H/(sum(H(:))+eps);
s1 = sum(H,1);
s2 = sum(H,2);
H  = H.*log2((H+eps)./(s2*s1+eps));
o = -sum(H(:));

global nrOfEvaluations
nrOfEvaluations = nrOfEvaluations+1;
return;





% ===========================================================================
% ===========================================================================
function display_results(targetInfo,objectInfo,x)
%
fig = spm_figure('FindWin','Graphics');
if isempty(fig), return; end;
set(0,'CurrentFigure',fig);
spm_figure('Clear','Graphics');

% Display text
%-----------------------------------------------------------------------
ax = axes('Position',[0.1 0.8 0.8 0.15],'Visible','off','Parent',fig);
text(0.5,0.7, 'Mutual Information Coregistration','FontSize',16,...
	'FontWeight','Bold','HorizontalAlignment','center','Parent',ax);

Q = inv(objectInfo.mat\spm_matrix(x(:)')*targetInfo.mat);
text(0,0.5, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),'Parent',ax);
text(0,0.3, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),'Parent',ax);
text(0,0.1, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),'Parent',ax);

% Display scatter-plots
%-----------------------------------------------------------------------
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
H   = spm_hist2(targetInfo.uint8,objectInfo.uint8,objectInfo.mat\targetInfo.mat,[1 1 1]);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
	'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
	'XTick',[],'YTick',[]);
title('Original Histogram','Parent',ax);
xlabel(spm_str_manip(targetInfo.fname,'k22'),'Parent',ax);
ylabel(spm_str_manip(objectInfo.fname,'k22'),'Parent',ax);

H   = spm_hist2(targetInfo.uint8,objectInfo.uint8,objectInfo.mat\spm_matrix(x(:)')*targetInfo.mat,[1 1 1]);
ax  = axes('Position',[0.6 0.5 0.35 0.3],'Visible','off','Parent',fig);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
	'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
	'XTick',[],'YTick',[]);
title('Final Histogram','Parent',ax);
xlabel(spm_str_manip(targetInfo.fname,'k22'),'Parent',ax);
ylabel(spm_str_manip(objectInfo.fname,'k22'),'Parent',ax);

% Display ortho-views
%-----------------------------------------------------------------------
spm_orthviews('Reset');
h1 = spm_orthviews('Image',targetInfo.fname,[0.01 0.01 .48 .49]);
h2 = spm_orthviews('Image',objectInfo.fname,[.51 0.01 .48 .49]);
global st
st.vols{h2}.premul = inv(spm_matrix(x(:)'));
spm_orthviews('Space',h1);

return;





% ===========================================================================
% ===========================================================================
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 & V.pinfo(1) == 2,
	mx = 255*V.pinfo(1) + V.pinfo(2);
	mn = V.pinfo(2);
else,
	fprintf('   Calculating maximum and minimum voxel values      ');
	mx = -Inf; mn =  Inf;
	for p=1:V.dim(3),
		img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
		mx  = max([max(img(:)) mx]);
		mn  = min([min(img(:)) mn]);
		fprintf('\b\b\b\b\b');
		fprintf('%-3d %%', round(100*p/V.dim(3)));
	end;
	fprintf('\n');
end;

fprintf('   Loading      ');
udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
for p=1:V.dim(3),
	img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
	udat(:,:,p) = uint8(round((img-mn)*((256-1)/(mx-mn))));
	fprintf('\b\b\b\b\b');
	fprintf('%-3d %%', round(100*p/V.dim(3)));
end;
fprintf('\n');
	
return;





% ===========================================================================
% ===========================================================================
function prettyEcho(objectName, targetName, others, free, startingParams)
%
fprintf('\n\n');
disp('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
fprintf('\n\n');
disp('                                        _                            ')
disp('                                       |_|                           ')
disp('           ___  __  __  ___     __  __  _  ___   ___  ____           ')
disp('          | __)|  \/  |/ __)   |  \/  || ||   ) | __)/  __)          ')
disp('          | __||      |\__ \   |      || ||   \ | __|| (_ |          ')
disp('          |___)|_|\/|_|(___/___|_|\/|_||_||_|\_\|___)\____|          ') 
fprintf('\n\n');
disp('               Koen Van Leemput       August 23, 2001                ');
fprintf('\n');
disp('                      koen.vanleemput@hut.fi                         ');
fprintf('\n');
disp('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
fprintf('\n\n');


fprintf('      objectName = %s\n', objectName);
fprintf('      targetName = %s\n', targetName);

fprintf('      others = [\n');
for i=1:size(others,1)
  fprintf('        ''%s'';\n', others(i,:));
end
fprintf('                  ]\n');

fprintf('      free = [')
for i=1:11
  fprintf('%d ', free(i))
end
fprintf('%d]\n', free(12))

fprintf('      startingParams = [')
for i=1:11
  fprintf('%f ', startingParams(i))
end
fprintf('%f]\n', startingParams(12))


fprintf('\n\n');
return





% ===========================================================================
% ===========================================================================
function printInterval(left, right, middle, middleValue)
%

fprintf('   [%15.10f   %15.10f  ]   %15.10f   %15.10f\n', ...
	left, right, middle, middleValue)
return





% ===========================================================================
% ===========================================================================
function printPosition(x)
%

fprintf('   translation (mm): %15.10f%15.10f%15.10f\n', x(1:3));
fprintf('   rotation (rad)  : %15.10f%15.10f%15.10f\n', x(4:6));
fprintf('   scaling         : %15.10f%15.10f%15.10f\n', x(7:9));
fprintf('   affine          : %15.10f%15.10f%15.10f\n', x(10:12));
fprintf('\n');
return





% ===========================================================================
% ===========================================================================
function uiForEms_mireg()
%
SPMid = spm('FnBanner',mfilename,'2.9');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','ems_mireg');

nsubjects = spm_input('number of subjects',1, 'e', 1);
if (nsubjects < 1)
  spm_figure('Clear','Interactive');
  return;
end

atlasNormalization = spm_input('Normalize to atlas?', 2, 'y/n', ...
    [1; 0]);

if atlasNormalization
  startingParams = [zeros(1,6) ones(1,3) zeros(1,3)];
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
  
  napList = [2 3 5 6 8 9 11 12];
  defaultNapNr = 4;
  nap = spm_input('# Affine Params?',3,'m',[
    ' 2 params (X & Y Translations)|'...
        ' 3 params (Translations)|'...
        ' 5 params (Translations + Pitch & Roll)|'...
        ' 6 params (Rigid Body)|'...
        ' 8 params (Rigid Body + X & Y Zoom)|'...
        ' 9 params (Rigid Body + Zooms)|'...
        '11 params (Rigid Body,  Zooms + X & Y Affine)|'...
        '12 params (Rigid Body,  Zooms & Affine)'
    ], napList, defaultNapNr);
  free = [ones(1,nap) zeros(1,12-nap)];
  
  startingParams = [0 0 0 0 0 0 1 1 1 0 0 0];
end


spm('Pointer','Watch');
for i=1:nsubjects
  set(spm_figure('FindWin','Interactive'),...
      'Name',['Coregistering subject ' num2str(i)],'Pointer','Watch');
  drawnow;
  eval(['targetName    =    targetName' num2str(i) ';']);
  eval(['objectName    =    objectName' num2str(i) ';']);
  eval(['others = others' num2str(i) ';']);
  
  ems_mireg(objectName, targetName, others, free, startingParams)
  
end
spm_figure('Clear',spm_figure('FindWin','Interactive'));
spm('FigName','ems_mireg: done',Finter,CmdLine);
spm('Pointer');
return;




