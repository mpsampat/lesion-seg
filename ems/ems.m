function R1=ems(Action)
% EMS: Automated segmentation of brain MR images
%_______________________________________________________________________
%    ___  __  __  ___ 
%   / __)(  \/  )/ __)  Expectation-Maximization Segmentation
%   ) _)  )    ( \__ \  Koen Van Leemput
%   \___)(_/\/\_)(___/  August 23, 2001
%_______________________________________________________________________
%
% 
%
% EMS is a freely available suite of MatLab functions and subroutines
% (with some externally compiled C routines) for fully-automated
% multi-spectral classification of brain tissues in Magnetic Resonance
% (MR) images. It is an add-on to the freely distributed SPM-package
% (Statistical Parametric Mapping, the Wellcome Department of
% Cognitive Neurology), and implements the methods described in the
% following papers: 
% 
% - K. Van Leemput, F. Maes, D. Vandermeulen, A. Colchester,
%   P. Suetens, Automated segmentation of multiple sclerosis lesions
%   by model outlier detection , IEEE transactions on medical imaging,
%   vol. 20, no. 8, pp. 677-688, August 2001
%
% - K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens, Automated
%   model-based tissue classification of MR images of the brain , IEEE
%   transactions on medical imaging, vol. 18, no. 10, pp. 897-908,
%   October 1999
%
% - K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens, Automated
%   model-based bias field correction of MR images of the brain , IEEE
%   transactions on medical imaging, vol. 18, no. 10, pp. 885-896,
%   October 1999
%
% - F. Maes, A. Collignon, D. Vandermeulen, G. Marchal, P. Suetens,
%   Multimodality image registration by maximization of mutual
%   information , IEEE transactions on Medical Imaging, vol. 16, no. 2,
%   pp. 187-198, April 1997 (in December 2000 recognized by the IEEE as
%   the journal's most frequently cited paper published in 1997)
%
%                           ----------------
%
% The EMS software is distributed free of charge for research
% purposes. It is supplied as is; no formal support or maintenance is
% provided. Additional information may be found at the EMS web site:
%       http://bilbo.esat.kuleuven.ac.be/web-pages/downloads/ems/
%                           ----------------
%
%              Please report bugs to <koen.vanleemput@hut.fi>
%
% ------------------------------------------------------------------------
% ems.m    Koen Van Leemput - August 23, 2001



emsVer    = 'EMS - August 23, 2001';

if nargin == 0, Action='Init'; end


if strcmp(lower(Action),lower('Init'))
%=======================================================================
ems_defaults;
global SPM_DIR 
addpath(SPM_DIR)
if isempty(spm_figure('FindWin','Menu'))
	spm('fmri')
	clc
else
	clc
end
ems('AsciiWelcome')
spm_help('!Disp','ems.m','','Graphics',ems('Ver'))
ems('CreateMenuWin')


elseif strcmp(lower(Action),lower('AsciiWelcome'))
%=======================================================================

disp( '    ___  __  __  ___                                         ')
disp( '   / __)(  \/  )/ __)  Expectation-Maximization Segmentation ')
disp( '   ) _)  )    ( \__ \  Koen Van Leemput                      ')
disp(['   \___)(_/\/\_)(___/  Version: ', ems('Ver')])



fprintf('\n')


elseif strcmp(lower(Action),lower('Ver'))
%=======================================================================
% ems('Ver')
R1 = emsVer;


elseif strcmp(lower(Action),lower('CreateMenuWin'))
%=======================================================================
close(findobj(get(0,'Children'),'Tag','EMS Menu'))

%-Open EMS menu window
%-----------------------------------------------------------------------

WS   = spm('WinScale');	
FS   = spm('FontSizes');
PF   = spm_platform('fonts');	
Rect = [508-200  429-150 400 310] .* WS;

F = figure('IntegerHandle','off',...
	'Name',sprintf('%s%s',ems('ver'),spm('GetUser',' (%s)')),...
	'NumberTitle','off',...
	'Tag','EMS Menu',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'MenuBar','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(12),...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(12),...
	'DefaultUicontrolInterruptible','on',...
	'Renderer','zbuffer',...
	'Visible','off');





%-Frames and text
%-----------------------------------------------------------------------
axes('Position',[0 0 100/400 280/300],'Visible','Off');

text(0.5,0.475,'EMS', 'FontName',spm_platform('Font','times'), ...
	 'FontSize',spm('FontSize', 72), 'FontWeight','Bold', ...
	 'Rotation',90, 'VerticalAlignment','middle', ...
	 'HorizontalAlignment','center', 'Color',[1 1 1]*.6);


axes('Position',[0 250/300 400/400 50/300],'Visible','Off');

text(0.5,0.5,'Expectation-Maximization Segmentation', ...
	  'FontName', spm_platform('Font','times'), ...
	  'FontSize', spm('FontSize', 17), ...
	  'FontAngle','Italic', 'FontWeight','Bold', ...
	  'VerticalAlignment','middle', ...
	  'HorizontalAlignment', 'center', 'Color',[1 1 1]*.6);

uicontrol(F,'Style','Frame','Position',[095 010 300 250].*WS,...
	'BackgroundColor',ems('Colour'));
uicontrol(F,'Style','Frame','Position',[105 060 135 190].*WS);
uicontrol(F,'Style','Frame','Position',[250 060 135 190].*WS);
uicontrol(F,'Style','Frame','Position',[105 020 280 30].*WS);
uicontrol(F,'Style','Text',...
	'String','Analysis',...
	'Position',[115 210 115 030].*WS,...
	'ForegroundColor','w','FontName',PF.times,'FontAngle','Italic');
uicontrol(F,'Style','Text',...
	'String','Tools',...
	'Position',[260 210 115 030].*WS,...
	'ForegroundColor','w','FontName',PF.times,'FontAngle','Italic');



%-Buttons to launch EMS functions
%-----------------------------------------------------------------------

global MIRIT_CMD
if exist(deblank(MIRIT_CMD))==2  
  uicontrol(F,'String','Register',...
	    'Position',[115 175 80 30].*WS,...
	    'ToolTipString','Mutual information based affine registration',...
	    'CallBack','ems_MIRIT');
else
  uicontrol(F,'String','Register',...
	    'Position',[115 175 80 30].*WS,...
	    'ToolTipString','Mutual information based affine registration',...
	    'CallBack','ems_mireg');
end


uicontrol(F,'String','?',...
	'Position',[205 175 25 30].*WS,...
	'CallBack','spm_help(''ems_MIRIT.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Segment',...
	'Position',[115 130 80 30].*WS,...
	'ToolTipString', ...
	       'Automated segmentation of normal brain MR images',...
	'CallBack','ems_segment');

uicontrol(F,'String','?',...
	'Position',[205 130 25 30].*WS,...
	'CallBack','spm_help(''ems_segment.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Lesions',...
	'Position',[115 085 80 30].*WS,...
	'ToolTipString', ...
	       'Automated segmentation of signal abnormalities in pathological brain MR images',...
	'CallBack','ems_lesions');

uicontrol(F,'String','?',...
	'Position',[205 085 25 30].*WS,...
	'CallBack','spm_help(''ems_lesions.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');




uicontrol(F,'String','Check atlas',...
	'Position',[260 195 80 20].*WS,...
	'ToolTipString','Check spatial position of images with respect to atlas',...
	'CallBack','ems_compareToAtlas');

uicontrol(F,'String','?',...
	'Position',[350 195 25 20].*WS,...
	'CallBack','spm_help(''ems_compareToAtlas.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Correct',...
	'Position',[260 165 80 20].*WS,...
	'ToolTipString','Correct images for estimated bias fields',...
	'CallBack','ems_getBiasCorrectedImage');

uicontrol(F,'String','?',...
	'Position',[350 165 25 20].*WS,...
	'CallBack','spm_help(''ems_getBiasCorrectedImage.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Resample',...
	'Position',[260 135 80 20].*WS,...
	'ToolTipString','Resample images',...
	'CallBack','ems_reslice');

uicontrol(F,'String','?',...
	'Position',[350 135 25 20].*WS,...
	'CallBack','spm_help(''ems_reslice.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Hard',...
	'Position',[260 105 80 20].*WS,...
	'ToolTipString','Assign voxels entirely to most probable class',...
	'CallBack','ems_makeHardSegmentation');

uicontrol(F,'String','?',...
	'Position',[350 105 25 20].*WS,...
	'CallBack','spm_help(''ems_makeHardSegmentation.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

u = uicontrol(F,'String','gipl2img',...
	'Position',[260 075 80 20].*WS,...
	'ToolTipString','Convert images from GIPL-format into SPM-format',...
	'CallBack','ems_gipl2img');

uHelp = uicontrol(F,'String','?',...
	'Position',[350 075 25 20].*WS,...
	'CallBack','spm_help(''ems_gipl2img.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

if exist('ems_gipl2img')~=2  
  set(u, 'enable', 'off')
  set(uHelp, 'enable', 'off')
end



uicontrol(F,'String','About EMS',...
	'Position',[115 025 80 020].*WS,...
	'ToolTipString','EMS man pages',...
	'CallBack','spm_help(''ems.m'')',...
	'ForegroundColor','g');

uicontrol(F,'String','EMSweb',...
	'Position',[205 025 80 020].*WS,...
	'ToolTipString', ...
	  'Launch web browser - http://bilbo.esat.kuleuven.ac.be/web-pages/downloads/ems/',...
        'CallBack','web(''http://bilbo.esat.kuleuven.ac.be/web-pages/downloads/ems/'');');

uicontrol(F,'String','Close EMS',...
	'Position',[295 025 80 020].*WS,...
	'ToolTipString','Close this window',...
        'CallBack','close(gcf)',...
	'ForegroundColor','r');


set(F,'Pointer','Arrow','Visible','on')




elseif strcmp(lower(Action),lower('Colour'))
%=======================================================================
% ems('Colour')
%-----------------------------------------------------------------------
R1 = [0.7,1.0,0.7];
%R1 = [0.8 0.8 1.0];



else
%=======================================================================
error('Unknown action string')

%=======================================================================
end