function varargout = AMRtestbed(varargin)
% AMRTESTBED MATLAB code for AMRtestbed.fig
% 
%      H = AMRTESTBED returns the handle to a new AMRTESTBED or the handle to
%      the existing singleton*.
%
%      AMRTESTBED('CALLBACK',hObject,eventData,handles,...) calls the
%      local function named CALLBACK in AMRTESTBED.M with the given input arguments.
%
%      AMRTESTBED('Property','Value',...) creates a new AMRTESTBED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AMRtestbed_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AMRtestbed_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AMRtestbed

% Last Modified by GUIDE v2.5 22-Dec-2018 21:40:41

%% List of callback function
% handles = pushbutton_inputSignal_Callback(hObject, eventdata, handles);    
% pushbutton_viewWaveform_Callback(hObject, eventdata, handles)
% pushbutton_viewConstellation_Callback(hObject, eventdata, handles)
% handles = pushbutton_viewSpectrum_Callback(hObject, eventdata, handles);    
% handles = pushbutton_addNoise_Callback(hObject, eventdata, handles);
% function pushbutton_applyFading_Callback(hObject, ~, handles)
% handles = pushbutton_coarseBandWidthEstimation_Callback(hObject, eventdata, handles);
% handles = pushbutton_coarseSNRestimation_Callback(hObject, eventdata, handles);
% handles = pushbutton_coarsefreqOffsetEstimation_Callback(hObject, eventdata, handles);
% handles = pushbutton_lowPassFilterling_Callback(hObject, eventdata, handles);
% handles = pushbutton_coarseSymRateEstimation_Callback(hObject, eventdata, handles);
% handles = pushbutton_reSampling_Callback(hObject, eventdata, handles);
% handles = pushbutton_coarseToneSpacing_Callback(hObject, eventdata, handles);
% handles = pushbutton_featureExtraction_Callback(hObject, eventdata, handles);
% handles = pushbutton_fusionCenter_Callback(hObject, eventdata, handles);

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AMRtestbed_OpeningFcn, ...
                   'gui_OutputFcn',  @AMRtestbed_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before AMRtestbed is made visible.
function handles = AMRtestbed_OpeningFcn(hObject, ~, handles, varargin)

%% Set default parameters
handles.opMode             = 0;  % operation mode : auto(1), manual(0)
handles.previousDecimationFactor = 1024;  
handles.sam_freq           = 140e6/handles.previousDecimationFactor;
handles.sam_period         = 1 / handles.sam_freq;
handles.rF                 = 0.35;    % roll-off factor
handles.FFTlength          = 2^13;    % 8192, The length of periodgram = FFT length, for Bandwidth
handles.samplesPerSegment  = 2^14;

handles.SNR_dB             = 15;
handles.limitation         = 120;          % Unused
handles.displayRange       = 1:1500;
handles.viewStartIdx       = 1;
handles.viewEndIdx         = 1500;
handles.traceType          = 'ClearWrite';  
handles.xdB                = 21;            % Unused, heuristic value
handles.oversampleFactor   = 26;            % for cyclic cumulant
handles.refDataType        = 'E4438C';
handles.currentPlot        = 0;
handles.freqResolution     = handles.sam_freq / handles.FFTlength;
% For classification between linear and non-linear modulation signal,
% refer to testLinearVSnonLinear.m
handles.sigmaAThreshold    = 0.26;  % edited by AWH at 10-Sep-2015        
handles.lowerConstantBW    = 0.5; 
% usually the estimated BW is greater than symbol rate due to the effect of rc filter
% thus to reduce the upper bound is recommended.
handles.upperConstantBW    = 1.5; 
handles.cutOffExpansionRatio = 1.3; % lowpass filtering bandwith constant
handles.MINIMUM_HIST_COUNT  = 13; % the counts of histogram in bandwidth estimation

handles.cutOffFreqWideBandLPF                = 100e3;  % Unused
handles.minimumSymbolRate                    = 900;
handles.maximumSymbolRate                    = 45e3;
handles.maxixumBandwidthForLinearModulation  = 45e3 * (1+0.5); % heuristic value
handles.AM_SIGNALFLOOR_THRESHOLD             = 10;
% handles.SUM_SIGNALFLOOR_UPPER_REGION         = 20; % Unused
% handles.RATIO_THRESHOLD                      = 5;  % Unused
% handles.SYMBOL_RATE_MAG_SPECTRUM_THRESHOLD   = 5;  % Unused

handles.nthOrder      = 2;
handles.qConj         = 0;
handles.dBFlag        = 1;
handles.spsLowerBound = 18;
handles.spsUpperBound = 35;
handles.bandwidthConstant = 0.7; % 11-Sep-2015
% Effect on fading channel
handles.maxDopplerShift = 107;
handles.fadingCarierFreq = 900;
handles.averageSpeedOfObject = 128.74; % km/h
handles.fadingType = 'Rayleigh';

 handles.isDehoppingSig = 0;
%% Remind that these parameter should be copy to reset callback

% set parameter
% markerSet={'+','o','*','.','x','s','d','^','v','>','<','p','h'};
% colorSet={'','g','b','c','m','k','r','g','b','c','m','k','r','g','b'};
handles.Lmod       = {'AM','FM','2ASK','2FSK', '4FSK', '8FSK','BPSK', 'QPSK',...
                          '8PSK','QAM'};
handles.numOfLienarModScheme   = {'QPSK','8PSK','32QAM'}; % except BPSK, 2ASK                         
% handles.Lmod       = {'AM','FM','2ASK','2FSK', '4FSK', '8FSK','BPSK', 'QPSK',...
%                           '8PSK','16QAM','32QAM','64QAM' };
% handles.Nmod         = {'doubleAm' ,'fm','ask2','fsk2','fsk4','fsk8','bpsk', ...
%                          'qpsk','psk8','qam16','qam32','qam64'};
handles.lineSet      = {'.r','.g','.b','.c','.m','.k','.r','.g','.b','.c','.m','.k','.r','.g'};
handles.markerSet      = {'o','+','*','.','s','d','^','v','>','p','h','<','o','+'};

% handles.Fset =  {'\gamma_{max}','\sigma_{aa}', '\sigma_{af}','\sigma_{a}' ,...
%                 '|C61|',...
%                 '|CC40|','|CC61|','|CC80|','|CC82|',...
%                 '|\mu_2|','Variance- \sigma^2'};
% handles.Fset =  {'\sigma_{aa}', '\sigma_{af}','\sigma_{a}' ,... // for KNN
%                 '|CC40|','|CC61|','|CC80|','|CC82|'};
handles.Fset =  {'\sigma_{a}','|CC40|','|CC20|'};  % for KNN
                
handles.refModIdx       = 1 : length(handles.numOfLienarModScheme);
handles.numOfFeature    = length(handles.Fset);

set(handles.edit_rF,            'String',num2str(handles.rF) );
set(handles.edit_samFreq,       'String',num2str(handles.sam_freq));
set(handles.edit_dispSNR,       'String','40');
set(handles.edit_viewStartIdx,  'String',num2str(handles.viewStartIdx));
set(handles.edit_viewEndIdx,    'String',num2str(handles.viewEndIdx));
set(handles.edit_freqResolution,'String',num2str(handles.freqResolution));
set(handles.edit_dispSNR,       'String',handles.SNR_dB);
set(handles.edit_xdB,           'String',num2str(handles.xdB));
set(handles.edit_samplesPerSegment,'String',num2str(handles.samplesPerSegment));
set(handles.popupmenu_decimationFactor,'Value',nextpow2(handles.previousDecimationFactor)-4);
set(handles.popupmenu_spectrumTrace,'Value',1);
set(handles.popupmenu_fftLength,'Value',5);    % 1024
set(handles.popupmenu_spectrumTrace,'Value',2);
set(handles.edit_bandwidthConstant,'String',num2str(handles.bandwidthConstant));
set(handles.checkbox_dBFlag,'Value',1);

set(handles.edit_fadingCarrierFreq,'String',num2str(handles.fadingCarierFreq));
set(handles.edit_maxDopplerShift,'String',num2str(handles.maxDopplerShift));
set(handles.edit_fadingObjSpeed,'String',num2str(handles.averageSpeedOfObject));

set(handles.slider_snr,'Value',handles.SNR_dB);
set(handles.edit_ddcDecimationRate,'Enable','on');
%% disable functions

set(handles.edit_viewEndIdx,'Enable','off');
set(handles.edit_viewStartIdx,'Enable','off');

set(handles.pushbutton_viewSpectrum,'Enable','off');
set(handles.pushbutton_viewWaveform,'Enable','off');
set(handles.pushbutton_viewConstellation,'Enable','off');
set(handles.uipushtool_viewForwardFig,'Enable','off');
set(handles.uipushtool_viewBackFig,'Enable','off');

set(handles.pushbutton_inputSignal,'Enable','on');
set(handles.pushbutton_nthPower,'Enable','off');
set(handles.popupmenu_nthPower,'Enable','off');
set(handles.pushbutton_addNoise,'Enable','off');
set(handles.pushbutton_applyFading,'Enable','off');

set(handles.pushbutton_coarseBandWidthEstimation,'Enable','off');
set(handles.pushbutton_coarsefreqOffsetEstimation,'Enable','off');
set(handles.pushbutton_coarseSNRestimation,'Enable','off');
set(handles.pushbutton_coarseSymRateEstimation,'Enable','off');
set(handles.pushbutton_coarseToneSpacing,'Enable','off');

set(handles.pushbutton_reSampling,'Enable','off');
set(handles.slider_snr,'Enable','on');
set(handles.pushbutton_addNoise,'Enable','off');

set(handles.pushbutton_lowPassFilterling,'Enable','off');
set(handles.pushbutton_featureExtraction,'Enable','off');
set(handles.pushbutton_fusionCenter,     'Enable','off');

set(handles.pushbutton_saveData,'Enable','off');
set(handles.pushbutton_makeCHeader,'Enable','off');
set(handles.pushbutton_demod,'Enable','off');
%% Disable Default setting
set(handles.popupmenu_decimationFactor,'Enable','off');

% Choose default command line output for AMRtestbed
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function varargout = AMRtestbed_OutputFcn(hObject, ~, handles) 
varargout{1} = handles.output;
varargout{2} = handles;

function handles = pushbutton_restart_Callback(hObject, ~, handles)
close(gcbf);
AMRtestbed;

%%%%%%%%%%%%%%%%%%%%%%%%%%   >Default setting<   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_filename_CreateFcn(hObject, eventdata, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_samplesPerSegment_Callback(hObject, eventdata, handles)
samplesPerSegment = str2double( get(hObject,'String') ); 
handles.samplesPerSegment = 2^nextpow2(samplesPerSegment);
set(hObject,'String',num2str(handles.samplesPerSegment));
guidata(hObject,handles);

function edit_samplesPerSegment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_decimationFactor_Callback(hObject, eventdata, handles)

contents                 = cellstr(get(hObject,'String'));
currentDecimationFactor  = str2double( contents{get(hObject,'Value')} );
previousDecimationFactor = handles.previousDecimationFactor;
rx_sig             = handles.rx_sig;
FFTlength                = handles.FFTlength;

power      = nextpow2(currentDecimationFactor);
integerNum = currentDecimationFactor / previousDecimationFactor;
%% update states
if integerNum == 1      % None
%     rx_sig = rx_sig;
elseif integerNum > 1   % Downsampling
%     rx_sig = decimate(basic_rx_sig,integerNum);
    rx_sig = downsample(rx_sig,integerNum);
    
elseif integerNum < 1   % Upsampling
%     rx_sig = interp(basic_rx_sig,1/integerNum);
    rx_sig = interp(rx_sig,1/integerNum);
end

sam_freq =  140e6/currentDecimationFactor;
sigLen   =  length(rx_sig);
%% update GUI states
cla(handles.axes_monitoring,'reset');
legend(handles.axes_monitoring,'off');

set(handles.edit_samFreq,'String',num2str(sam_freq));
set(handles.edit_inputNumSamples,'String',num2str(sigLen));
set(handles.edit_viewEndIdx,'String',num2str(sigLen));
% set(hObject,'Value',power-4);

%% update handles structure 
handles.sam_freq          = sam_freq;
handles.sam_period        = 1 / handles.sam_freq;
handles.rx_sig            = rx_sig;
handles.viewEndIdx        = sigLen;
handles.freqResolution    = sam_freq / FFTlength;    
handles.previousDecimationFactor  = currentDecimationFactor;
guidata(hObject,handles);

function edit_dispSNR_Callback(hObject, eventdata, handles)
handles.SNR_dB = str2double( get(hObject,'String') );
guidata(hObject,handles);

function edit_dispSNR_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','10');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_rF_Callback(hObject, eventdata, handles)
handles.rF = str2double( get(hObject,'String') );
guidata(hObject,handles);
function edit_rF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_samFreq_Callback(hObject, eventdata, handles)
handles.sam_freq = str2double( get(hObject,'String') );
guidata(hObject,handles);
function edit_samFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_inputNumSamples_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_inputNumSymbols_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewStartIdx_Callback(hObject, ~, handles)
handles.viewStartIdx = str2double( get(hObject,'String') );
if handles.viewStartIdx <= 1
    handles.viewStartIdx = 1;
end
set(hObject,'String', num2str(handles.viewStartIdx));
guidata(hObject,handles);
function edit_viewStartIdx_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_viewEndIdx_Callback(hObject, eventdata, handles)
handles.viewEndIdx = str2double( get(hObject,'String') );
sigMaxLength = length(handles.basic_rx_sig);
if handles.viewEndIdx >= sigMaxLength
    handles.viewEndIdx = sigMaxLength;
end
set(hObject,'String', num2str(handles.viewEndIdx));

guidata(hObject,handles);
function popupmenu_spectrumTrace_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String')); % returns popupmenu contents as cell array
traceType = (contents{get(hObject,'Value')}); % returns selected item from popupmenu

handles.traceType = traceType;
guidata(hObject,handles);
function popupmenu_spectrumTrace_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewEndIdx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_freqResolution_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_refDataType_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));  
refDataType = contents{get(hObject,'Value')};         
handles.refDataType = refDataType;
guidata(hObject,handles);
function popupmenu_refDataType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_decimationFactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_upsamplingRatio_Callback(hObject, eventdata, handles)

% function handles=popupmenu_decimationFactor_Callback(hObject, ~, handles)
% %
% % Perform interpolation or decimation 
% %
% decimationFactor = handles.decimationFactor;
% defaultDecimationFactor = handles.defaultDecimationFactor;
% basic_rx_sig = handles.basic_rx_sig;
% FFTlength = handles.FFTlength;
% 
% contents = cellstr(get(hObject,'String')); % returns popupmenu_decimationFactor contents as cell array
% selectedIntNumber = str2double(contents{get(hObject,'Value')}); % returns selected item from popupmenu_decimationFactor
% 
% currentDecimationFactor = selectedIntNumber*defaultDecimationFactor;
% 
% %% update states
% if currentDecimationFactor ~= decimationFactor    
%     rx_sig = decimate(basic_rx_sig,selectedIntNumber);
% 
%     set(handles.edit_samFreq,'String',num2str(  140e6 / currentDecimationFactor ));
%     set(handles.edit_inputNumSamples,'String',num2str(length(rx_sig)));
%     set(handles.edit_viewDecimationFactor,'String',num2str( currentDecimationFactor));
%     
%     sam_freq = 140e6 / currentDecimationFactor;
%     
%    
%    %% save data      
%     handles.decimationFactor = currentDecimationFactor;
%     handles.rx_sig = rx_sig;
%     handles.viewEndIdx = length(rx_sig);
%     
%     handles.sam_freq =  sam_freq;
%     handles.freqResolution = sam_freq / FFTlength;
% end
% guidata(hObject,handles);

function popupmenu_fftLength_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.FFTlength = str2double( contents{get(hObject,'Value')} );
% Set frequency resolution
set(handles.edit_freqResolution,'String',num2str(handles.sam_freq/handles.FFTlength));
guidata(hObject,handles);

function popupmenu_fftLength_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_xdB_Callback(hObject, ~, handles)
xdB = str2double(get(hObject,'String'));
handles.xdB = xdB;
guidata(hObject,handles);

function edit_xdB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_bandwidthConstant_Callback(hObject, eventdata, handles)
handles.bandwidthConstant = str2double(get(hObject,'String'));
guidata(hObject,handles);

function edit_bandwidthConstant_CreateFcn(hObject, eventdata, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_saveData_Callback(~, eventdata, handles)
filename = handles.filename;
rx_sig   = handles.rx_sig;
viewStartIdx  = handles.viewStartIdx;
viewEndIdx    = handles.viewEndIdx;

rx_sig = rx_sig( viewStartIdx : viewEndIdx);
%% save .mat file
% startIdx = strfind(filename,'.mat');
% filename(startIdx:end) = [];
% save(strcat(filename,'_conditioned','.mat'),'rx_sig');
startIdx = strfind(filename,'.diq');
filename(startIdx:end) = [];
save(strcat(filename,'_numIQ_',num2str(viewEndIdx-viewStartIdx+1), ...
    '_SNR_',num2str(handles.SNR_dB),'_date_',datestr(now,'yymmdd'),'.mat'),'rx_sig');

% function pushbutton_amrTraining_Callback(hObject, ~, handles)
% refDataType = handles.refDataType;
% 
% refData_generation(refDataType,handles)

function radiobutton_auto_Callback(hObject, eventdata, handles)
function radiobutton_manual_Callback(hObject, eventdata, handles)
function uipanel_opMode_SelectionChangeFcn(hObject, eventdata, handles)

switch hObject
    case handles.radiobutton_auto
        handles.opMode = 1;
    case handles.radiobutton_manual
        handles.opMode = 0;
end

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%   <Default setting>   %%%%%%%%%%%%%%%%%%%%%%%%
function handles = pushbutton_inputSignal_Callback(hObject, eventdata, handles)
updateLoggingMsg('While reading the files...',handles);
minIQinHop = 9000;
isDehoppingSig = 0;

cd(pwd);

%% Note
%  주 스펙트럼외에 부수적인 스펙트럼이 생긴 경우 IQ 데이터 파형을 살펴보고 
%  파형이 급격하게 변하는 부분에서 파형의 위 아래가 바뀌었는지 살펴볼 것.  
%  원인 : overflow
%  해결 : 신호 전력을 낮춘다. 


if isfield(handles,'sigName')
    filename = handles.sigName;
else
    [filename,filepath]=uigetfile({'*.*','All Files'},...
        'Select Data File 1');
end
%% read data
if filename ~= 0
    set(handles.edit_filename,'String',filename);
    
    [path,name,ext] = fileparts(filename);
    switch ext
        case {'.diq','.bin'} 
            fileID = fopen(filename,'r');
            if fileID > 0
                
                [initial_uint32_Data, count] = fread( fileID, 'uint32=>uint32');
                fclose(fileID);

                % restrict the number of data samples
%                 endNum = 2; % 20
%                 sigLenLimitStart    = 1;            
%                 sigLenLimitEnd      = floor(length( initial_uint32_Data ) / endNum);       
%                 initial_uint32_Data = initial_uint32_Data(sigLenLimitStart:sigLenLimitEnd);
                
                sigLenLimitStart    = 1;            
                sigLenLimitEnd      = length(initial_uint32_Data);       
                initial_uint32_Data = initial_uint32_Data(sigLenLimitStart:sigLenLimitEnd);
                
                % Discard first 16byte --> cancel, mis-understand
                uint32_Data = initial_uint32_Data(1:end);
                clear initial_uint32_Data;

                sigLen = sigLenLimitEnd; %length(uint32_Data ); 
                % Discard remaining byte
                uint32_Data = uint32_Data(1:(sigLen-mod(sigLen,4))); 
                % Flip data by 16byte 
                temp = reshape(uint32_Data,4,[]);
                uint32_Data = temp(4:-1:1,:);
                clear temp;    

                %%%%%%%%%%%%%%%%%%%%%%%% IQ data structure %%%%%%%%%%%%%%%%%%%%%%%%
                % Type      |  header  |        Inphase      |       Quadrature   |     
                % Bit index |  31   30 |     29     ~    15  |   14      ~      0 |           
                %           | TOA : 11 | (sign bit)          | (sign bit)         |
                %           | IQ  : 01 |                     |                    |
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
                
               %% Matlab style coding
                % Search TOA
                headerBits = bitand(uint32_Data,hex2dec('c0000000'));
%                 
                logicIdx = headerBits == hex2dec('c0000000');
                number = 1 : sigLen;
                TOAIdx = number(logicIdx);
                logic = diff(TOAIdx) >= minIQinHop; 
                Idx = number(logic);
                beginTOAidx = TOAIdx(Idx);
                endTOAidx = TOAIdx(Idx+1);
                numTOA = length(Idx);
%                 figure; stem(headerBits(:)); hold on; %xlim([0 3000]);
%                 len = length(beginTOAidx);
%                 plot(beginTOAidx,3221225472*ones(1,len),'r*'); %xlim([0 3000]);
%                 plot(endTOAidx,3221225472*ones(1,len),'go'); %xlim([0 3000]);
                               
%                 if ~isempty(Idx)
                    if sigLen - TOAIdx(end) >= minIQinHop;
                        beginTOAidx(1) = TOAIdx(end);   
                        endTOAidx(1) = sigLen;            
                    end
                     numTOA = 1;
%                 else
%                     disp('fail to detect the TOA');

%                 end

                offset = 10; rx_sig = 0; totniq = 0;
        
                for i = 1 : numTOA
                    tmp = uint32_Data(beginTOAidx(i)+offset : endTOAidx(i)-offset);
                    headerBits = bitand(tmp,hex2dec('c0000000'));
                    logicIdx = headerBits ~= hex2dec('40000000'); %                      
                    tmp(logicIdx) = [];

                    % Inphase
                    oasis = bitand( tmp,hex2dec('3fff8000'));
                    i_tmp = bitshift(oasis,-15); clear oasis;

                    % Quadrature
                    q_tmp = bitand( tmp,hex2dec('00007fff'));

                    % two's complement
                    maskingHex4000 = hex2dec('4000'); % sign bit
                    masking = hex2dec('3fff'); 
                    % find negative value
                    comp_I_Idx = bitand(i_tmp,maskingHex4000) == maskingHex4000;

                    inPhase = zeros(1,length(comp_I_Idx));

                    % Take two's complement
                    tmp = (masking - bitand( i_tmp(comp_I_Idx),masking )) +1; 
                    % Save value as double type in MATLAB
                    inPhase(comp_I_Idx) = -double(tmp);  
                    tmp = bitand(i_tmp(~comp_I_Idx),masking);
                    inPhase(~comp_I_Idx) = double(tmp);
                    clear comp_I_Idx;      

                    comp_Q_Idx = bitand(q_tmp,maskingHex4000) == maskingHex4000;

                    quadrature = zeros(1,length(comp_Q_Idx));
                    tmp = (masking - bitand( q_tmp(comp_Q_Idx),masking )) +1;
                    quadrature(comp_Q_Idx) =  -double(tmp);
                    tmp = bitand(q_tmp(~comp_Q_Idx),masking);
                    quadrature(~comp_Q_Idx) = double(tmp);
                    clear comp_Q_Idx;

                    amrObj(i).iq = inPhase + 1i*quadrature;
                    amrObj(i).niq = length(inPhase);
                    amrObj(i).IQbeginIdx = beginTOAidx(i)+offset;
                    amrObj(i).IQendIdx   = endTOAidx(i)-offset;

                    rx_sig = [ rx_sig amrObj(i).iq];
                    totniq = totniq + amrObj(i).niq;
                    clear inPhase; clear quadrature;
                end
                sigLenLimitEnd = length(rx_sig);
                     
                % Is the received signal FH?
                if numTOA < 10
                    isDehoppingSig = 0;
                    %normalization
                    rx_sig = sqrt(length(rx_sig)) * rx_sig / norm(rx_sig,2);  
                    updateLoggingMsg('Signal power is normalized to unit',handles)

                elseif numTOA > 10
                    isDehoppingSig = 1;
                end                
              
            else
                updateLoggingMsg('Move the exe file to the folder including diq data',handles);

            end
        case '.mat'
            sig = load(filename);
            fieldname = fieldnames(sig);
            dd =  strcat('sig.',fieldname);
            rx_sig = eval(dd{1});
                            
            sigLenLimit = length(rx_sig);
            sigLenLimitStart = 1;
            sigLenLimitEnd = sigLenLimit;
            isDehoppingSig = 0;
            handles.viewEndIdx = sigLenLimit;            
        case '.fig'
            % Open FIG-file with guide command.
            guide(filename)
            rx_sig = [1 2];
            isDehoppingSig = 0;
        case {'.wav', '.mp3'}
             % When size of a file is too large, 
             % restrict the number of loaded samples
  
             [rx_sig, handles.sam_freq] = ...
                  audioread(filename);             
            % Normalization
            rx_sig = sqrt(length(rx_sig)) * rx_sig / norm(rx_sig,2);  
            
            sigLenLimit = length(rx_sig);
            sigLenLimitStart = 1;
            sigLenLimitEnd = sigLenLimit;   
            isDehoppingSig = 0;
        otherwise
            try
                % Use open for other file types.
                open(filename)
            catch ex
                errordlg(...
                  ex.getReport('basic'),'File Type Error','modal')
            end

    end
    % Calculate FFT length
    avgNum = sigLenLimitEnd / 8;
    NFFT = 2.^nextpow2(avgNum);

    if abs(NFFT - avgNum) > abs(NFFT/2 - avgNum)
        NFFT = NFFT/2;
    end
    if NFFT > 8192
        NFFT = 8192;
    elseif NFFT < 128
        NFFT = 128;
    end
    
    nElementNum = nextpow2(NFFT) - 5; % mapping number
    set(handles.edit_inputNumSamples,'String',num2str(sigLenLimitEnd-sigLenLimitStart)); 
    set(handles.edit_viewEndIdx,'String',num2str(sigLenLimitEnd));
    set(handles.popupmenu_fftLength,'Value',nElementNum);
    updateLoggingMsg(filename,handles);
    
    %% save data
    handles.FFTlength = NFFT;
    handles.filename = filename;
    handles.rx_sig   = rx_sig;
    handles.basic_rx_sig = rx_sig;
    handles.viewEndIdx = sigLenLimitEnd;
    handles.isDehoppingSig = isDehoppingSig;
    guidata(hObject, handles);
    
%     rx_sig = rx_sig(1:85000);
%     startIdx = strfind(filename,'.diq');
%     filename(startIdx:end) = [];
%     filename = [filename,'.mat'];
%     save(filename,'rx_sig');
%     
%     makecheader_141118(filename);
    %% Auto mode
    if handles.opMode
          handles = pushbutton_addNoise_Callback(hObject, eventdata, handles);
          handles = pushbutton_coarseBandWidthEstimation_Callback(hObject, eventdata, handles);
          handles = pushbutton_coarseSNRestimation_Callback(hObject, eventdata, handles);
          handles = pushbutton_coarsefreqOffsetEstimation_Callback(hObject, eventdata, handles);
          handles = pushbutton_lowPassFilterling_Callback(hObject, eventdata, handles);
          handles = pushbutton_coarseSymRateEstimation_Callback(hObject, eventdata, handles);
          handles = pushbutton_coarseToneSpacing_Callback(hObject, eventdata, handles);
          handles = pushbutton_reSampling_Callback(hObject, eventdata, handles);
          handles = pushbutton_featureExtraction_Callback(hObject, eventdata, handles);
          handles = pushbutton_fusionCenter_Callback(hObject, eventdata, handles);
    else
    
        %% Update GUI states
        set(handles.popupmenu_decimationFactor,'Enable','on');

        set(handles.edit_viewEndIdx,'Enable','on');
        set(handles.edit_viewStartIdx,'Enable','on');

        set(handles.pushbutton_viewSpectrum,'Enable','on');
        set(handles.pushbutton_viewWaveform,'Enable','on');
        set(handles.pushbutton_viewConstellation,'Enable','on');

        set(handles.pushbutton_coarseBandWidthEstimation,'Enable','on');
        set(handles.pushbutton_nthPower,'Enable','on');
        set(handles.popupmenu_nthPower,'Enable','on');    
        
        set(handles.pushbutton_saveData,'Enable','on');
        set(handles.pushbutton_makeCHeader,'Enable','on');

        set(handles.pushbutton_inputSignal,'BackgroundColor',[222 235 250]/255);

       %% view spectrum
        handles=pushbutton_viewSpectrum_Callback(handles.pushbutton_viewSpectrum, eventdata, handles);
%         envelope=abs(rx_sig);                   % find the envelope 
%         myAmp = envelope-mean(envelope);
%         
%         nn = 3000;
%         xx = (0: nn-1); %/sam_freq;
%         xx2 = xx(1:end-1);
%         
%         myAmp = myAmp(1:nn);
%         myPhase = unwrap(angle(rx_sig(1:nn)));      
%         myFreq = diff(myPhase); 
%         myFreq = smooth(myFreq,20,'moving');  
        
        %% Requist for SHS
%         figure('name','Intantaneous information'); 
%         set(gcf,'Position',[500. 200.0 1065, 609]);
%         subplot(4,1,1), plot(xx,real(rx_sig(1:nn)),'-r.',xx,imag(rx_sig(1:nn)));title('IQ plot'); grid on;
%         subplot(4,1,2), plot(xx,myAmp,'-b'); xlabel('Samples'); ylabel('i amplitude');grid on;
%         subplot(4,1,3), plot(xx,myPhase/10,'-b'); xlabel('Samples'); ylabel('i phase');grid on;
%         subplot(4,1,4), plot(xx2,(myFreq),'-r.'); xlabel('Samples'); ylabel('i frequency');grid on;
       %% plot IQdata   
%     cla(handles.axes_monitoring,'reset');
%     legend(handles.axes_monitoring,'off');
%     hold(handles.axes_monitoring, 'on');
% 
%     if length(displayRange) > length(rx_sig)
%         displayRange = 1 : length(rx_sig);
%     end
% 
%     tt = (0 : length(displayRange)-1) * handles.sam_period; 
%     
%     plot(tt,real(rx_sig(displayRange)),'-r','parent',handles.axes_monitoring);
%     if ~isreal(rx_sig)
%         plot(tt,imag(rx_sig(displayRange)),'-b','parent',handles.axes_monitoring);
%         legend(handles.axes_monitoring,{'real','imag'})
%     end
%     xlabel(handles.axes_monitoring,'time[s]');
%     grid(handles.axes_monitoring,'on');
    
    end
end


function pushbutton_viewConstellation_Callback(hObject, eventdata, handles)
updateLoggingMsg('While processing constellation plot, please wait ...',handles);    
try
    
%     if isfield(handles,'noisy_rx_sig')
%         rx_sig = handles.noisy_rx_sig;
%     else
        rx_sig = handles.rx_sig;
%     end
    
    viewEndIdx = handles.viewEndIdx;
    viewStartIdx = handles.viewStartIdx;

    displayRange = viewStartIdx : viewEndIdx; 

    if ~handles.opMode
     %% plot IQdata
        cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off');
        hold(handles.axes_monitoring, 'on');

        plot(real(rx_sig(displayRange)),imag(rx_sig(displayRange)),'.','parent',handles.axes_monitoring);

        xlabel(handles.axes_monitoring,'In-Phase(real part)');
        ylabel(handles.axes_monitoring,'Quadrature(Imaginary part)');

        grid(handles.axes_monitoring,'on');

        updateLoggingMsg('Complete constellation plot',handles); 
    end
catch err
    updateLoggingMsg('Error in pushbutton_viewConstellation_Callback, Select the input signal',handles);
    rethrow(err);    
end

function pushbutton_viewWaveform_Callback(hObject, eventdata, handles)
updateLoggingMsg('While processing waveform plot, please wait ...',handles);    
try

%     if isfield(handles,'noisy_rx_sig')
%         rx_sig = handles.noisy_rx_sig;
%     else
        rx_sig = handles.rx_sig;
%     end

    viewEndIdx = handles.viewEndIdx;
    viewStartIdx = handles.viewStartIdx;

%     snr_dB = floor( get(hObject,'Value') );
%     set(handles.edit_dispSNR,'String',num2str(snr_dB));
    displayRange = viewStartIdx : viewEndIdx; %length(rx_sig); %handles.displayRange;

    if ~handles.opMode
        %% plot IQdata
        cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off');
        hold(handles.axes_monitoring, 'on');

        tt = (0 : length(displayRange)-1) * handles.sam_period; 
        plot(tt,real(rx_sig(displayRange)),'-r','parent',handles.axes_monitoring);
        if ~isreal(rx_sig(displayRange))
            plot(tt,imag(rx_sig(displayRange)),'-b','parent',handles.axes_monitoring);
            legend(handles.axes_monitoring,{'real','imag'})
        end
        xlabel(handles.axes_monitoring,'time[s]');
        grid(handles.axes_monitoring,'on');
        updateLoggingMsg('Complete waveform plot, please wait ...',handles);    
    end
catch err    
    updateLoggingMsg('Error in pushbutton_viewWaveform_Callback',handles);
    rethrow(err);
end

function checkbox_dBFlag_Callback(hObject, eventdata, handles)
handles.dBFlag = get(hObject,'Value');
guidata(hObject,handles);

function handles=pushbutton_viewSpectrum_Callback(hObject, eventdata, handles)
updateLoggingMsg('While processing spectrum plot, please wait ...',handles);    
try
    rx_sig       = handles.rx_sig;

    FFTlength    = handles.FFTlength;             
    sam_freq     = handles.sam_freq;
    viewEndIdx   = handles.viewEndIdx;
    viewStartIdx = handles.viewStartIdx;
    traceType    = handles.traceType;
    dBFlag       = handles.dBFlag;
    
    rx_sig = rx_sig(viewStartIdx : viewEndIdx );
    sigLen = length(rx_sig);
    rx_sig = rx_sig(:);
  
    %% zero-padding
    if sigLen < FFTlength
        rx_sig = [rx_sig; zeros(FFTlength-(viewEndIdx-viewStartIdx+1) ,1) ];
    end

    %% Calculated power spectral density
    % [Pxx,~] = pwelch(rx_sig,window,noverlap,FFTlength, ...
    %     sam_freq,'twosided');
    window           = hanning( FFTlength );            % Time domain data tapering window
    noverlap         = FFTlength/2;                     % The number of samples for overlap
    freqResolution   = sam_freq / FFTlength;

    numberOfSegments = (length(rx_sig)-noverlap)./(FFTlength-noverlap);
    numberOfSegments = fix(numberOfSegments);
    
    LminusOverlap = FFTlength-noverlap;
    xStart = 1:LminusOverlap:numberOfSegments*LminusOverlap;
    xEnd   = xStart+FFTlength-1;

    % set frequency axis
    freqAxis = ( -FFTlength/2 : (FFTlength -1)/2 ) / FFTlength * sam_freq;
    
    Sxx = 0; Sxxk=0;
    accumulatedSxx=  zeros(numberOfSegments,FFTlength);
    updateLoggingMsg(['Trace type is ',traceType],handles); 
    
    if(numberOfSegments > 19) numberOfSegments = 19; end  %to increase the AMR speed, recommend odd number at 150610
        
    switch traceType
        case 'MaxHold'
            for i = 1:numberOfSegments     
                xp = rx_sig(xStart(i):xEnd(i));
                
                xw = xp.*window;              % Window the data                                         
                U = window'*window;           % compensates for the power of the window.
                Xx = fft(xw,FFTlength);
                Sxxk = Xx.*conj(Xx)/U;        % Auto spectrum.
                accumulatedSxx(i,:)  = Sxxk;  % MaxHold Trace in ITU-R SM.443-4
                % Too slow
%                 updateLoggingMsg(['Processing ',num2str(i),' segments'],handles);    
%                 plot(freqAxis,fftshift( 10*log10(Sxxk) ),'--','color',rand(1,3),'parent',handles.axes_monitoring);
%                 xlabel('Frequency (Hz)','parent',handles.axes_monitoring);
%                 ylabel('Magnitude (dB/Hz)','parent',handles.axes_monitoring);
%                 xlim(handles.axes_monitoring,[freqAxis(1), freqAxis(end)]);
%                 grid(handles.axes_monitoring,'on'); 
%                 hold(handles.axes_monitoring,'on')

            end     
            if size(accumulatedSxx,1) > 1
                maxHoldSxx = max(accumulatedSxx);        
            else
                maxHoldSxx = (accumulatedSxx);
            end        
            Pxx = maxHoldSxx;

        case 'ClearWrite'
            for i = 1:numberOfSegments            
                xp = rx_sig(xStart(i):xEnd(i));
                xw = xp.*window;                % Window the data
                U = window'*window;             % compensates for the power of the window.
                Xx = fft(xw,FFTlength);
                Sxxk = Xx.*conj(Xx)/U;          % Auto spectrum.
                Sxx  = Sxx + Sxxk;
%                 updateLoggingMsg(['Processing ',num2str(i),' segments'],handles);
%                 plot(freqAxis,fftshift( 10*log10(Sxxk) ),'--','color',rand(1,3),'parent',handles.axes_monitoring);
%                 xlabel('Frequency (Hz)','parent',handles.axes_monitoring);
%                 ylabel('Magnitude (dB/Hz)','parent',handles.axes_monitoring);
%                 xlim(handles.axes_monitoring,[freqAxis(1), freqAxis(end)]);
%                 grid(handles.axes_monitoring,'on'); 
%                 hold(handles.axes_monitoring,'on')
            end     
            len =length( rx_sig(xEnd(end)+1 : end) );
            zz = FFTlength-len;
            residual = [rx_sig(xEnd(end)+1 : end); zeros(zz,1)];
            xw = residual.*window;                % Window the data
            U = window'*window;             % compensates for the power of the window.
            Xx = fft(xw,FFTlength);
            Sxxk = Xx.*conj(Xx)/U;          % Auto spectrum.
            Sxx  = Sxx + Sxxk;
                
            Pxx = Sxx ./ numberOfSegments; % Average the sum of the periodograms
        otherwise
    end

    Pxx = fftshift(Pxx);          Pxx_dB =( 10*log10(Pxx) );

    %% save data
    handles.Pxx = Pxx;
    handles.Pxx_dB = Pxx_dB;
    handles.freqAxis = freqAxis;
    handles.freqResolution = freqResolution;
    guidata(hObject,handles);

    %% For visualization
%     normFreqAxis = freqAxis./(2*max(freqAxis));
%     figure;  plot(normFreqAxis, Pxx_dB); grid off; 
%     axis([min(normFreqAxis),max(normFreqAxis),min(Pxx_dB),max(Pxx_dB)]);
%     xlabel('f/fs')
    
%     figure; hist(Pxx_dB, 35);
%     figure('name','Intantaneous information'); 
%     nn = 1200; xx = 1:nn;
%     plot(xx,real(rx_sig(xx)),'r',xx,imag(rx_sig(xx))); grid on;
%     xlabel('Samples')
%     axis( [min(xx),max(xx),min(real(rx_sig(xx))),max(real(rx_sig(xx)))] ) ;
    % plot graph    
    if ~handles.opMode
        legend(handles.axes_monitoring,'off');
        cla(handles.axes_monitoring,'reset');

        switch dBFlag
            case 0
                plot(freqAxis, Pxx,'parent',handles.axes_monitoring);
                ylim(handles.axes_monitoring,[min(Pxx), max(Pxx)]);
                ylabel('Magnitude','parent',handles.axes_monitoring);
                %     xlim(handles.axes_monitoring,[freqAxis(1), freqAxis(end)]);
                %     xlim(handles.axes_monitoring,[-60e3,60e3]);             

            case 1
                plot(freqAxis, Pxx_dB,'parent',handles.axes_monitoring);
                ylim(handles.axes_monitoring,[min(Pxx_dB), max(Pxx_dB)]);
                ylabel('Magnitude (dB)','parent',handles.axes_monitoring);
                %     xlim(handles.axes_monitoring,[freqAxis(1), freqAxis(end)]);
                %     xlim(handles.axes_monitoring,[-60e3,60e3]);           
            otherwise
        end    
        xlim(handles.axes_monitoring,[min(freqAxis), max(freqAxis)]);        
        xlabel('Frequency (Hz)','parent',handles.axes_monitoring);  
        grid(handles.axes_monitoring,'on');  
        %% Update GUI states    
        updateLoggingMsg('Complete plot processing',handles);  
        % set(handles.pushbutton_viewSpectrum,'BackgroundColor',[222 235 250]/255);
    end
catch err    
    updateLoggingMsg('Error in pushbutton_viewSpectrum_Callback',handles);
    rethrow(err);
end

% function pushbutton_preLPF_Callback(hObject, eventdata, handles)
% updateLoggingMsg('While low-pass filtering, please wait ...',handles);  
% 
% sam_freq  = handles.sam_freq;
% rx_sig    = handles.rx_sig;
% 
% cutOffFreqWideBandLPF = handles.cutOffFreqWideBandLPF;
% 
% filterOrder = 128;  
% Wn = (cutOffFreqWideBandLPF / (sam_freq ))  ;
% 
% lowPassFilterCoefficients = fir1(filterOrder,2*Wn);
% lowPassFilterCoefficients = lowPassFilterCoefficients(:);
% 
% %% filtering method using filter function in MATLAB
% rx_sig = filter(lowPassFilterCoefficients,1,rx_sig);
% 
% %% update GUI states
% set(handles.pushbutton_coarseBandWidthEstimation,'Enable','on');
% set(hObject,'BackgroundColor',[222 235 250]/255);
% updateLoggingMsg('Complete low-pass filtering',handles);  
% cla(handles.axes_monitoring,'reset');
% 
% %% save handles structure
% handles.rx_sig	= rx_sig;
% guidata(hObject,handles);


function popupmenu_nthPower_Callback(hObject, eventdata, handles)
%% get data from popupmenu
contents = cellstr(get(hObject,'String'));
nthOrderqConj = contents{get(hObject,'Value')};

slashIdx = strfind(nthOrderqConj,'/');

nthOrder = str2double(nthOrderqConj(1:slashIdx-1));
qConj = str2double(nthOrderqConj(slashIdx+1:end));

handles.nthOrder = nthOrder;
handles.qConj = qConj;
guidata(hObject,handles);

function popupmenu_nthPower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_nthPower_Callback(hObject, eventdata, handles)
nthOrder = handles.nthOrder;
qConj = handles.qConj;
rx_sig = handles.rx_sig;

poweredSig = rx_sig.^(nthOrder-qConj).* conj(rx_sig).^(qConj);

% switch coarseModType
%     case {'2FSK','4FSK','8FSK'}
%         % Obtain instantaenous frequency of continuous phase modulation(CPM) signal
%         phi_phase = angle(rx_sig);  % Between +- pi        
%         % make unwrapped angle
%         phi_NL_phase=unwrap(phi_phase); 
%         f_frequency=diff(phi_NL_phase);
%         sig = f_frequency; 
%         sig = sig - mean(sig);
%     otherwise % Linear modulation such as ASK,PSK,and QAM
%         % Step 2 : Estimate (2,1) cyclic moment       
%         sig = rx_sig;
% end
% %% Oerder and H.Meyr algorithm = Squaring method in MATLAB
% [normMag, freqResolution, NFFT] = ...
%     calculateCyclicMoment(sig,2,1,handles);  
% % To complement fftshift in calculateCyclicMoment function            
% normMag = fftshift(normMag);
% % Obtain the detail frequency vector 
% fineFreqAxis = (1:NFFT/2) /NFFT * sam_freq ;
% 
% normalizedMag_dB = 10*log10(normMag);
% 
% % must be fs >= 4*BW !!!, 
% % search symbol rate in range from minimum symbol rate to maximum
% % symbol rate
% count=0; partOfNormalizedMag_dB=[]; partOfFreqAxies=[];    
% % Move spectrum near to zero
% for i = 1 : length(fineFreqAxis)
%    if minimumSymbolRate <= fineFreqAxis(i) && fineFreqAxis(i) <= maximumSymbolRate
%        count = count + 1;
%        partOfNormalizedMag_dB(count) = normalizedMag_dB(i);
%        partOfFreqAxies(count) = fineFreqAxis(i);       
%    end    
% end
% 
% % use a moving average filter with a Nos span to smooth all the data at once
% Nos = round( NFFT * 0.1 ); %filter span
% spectralTrend = smooth(partOfNormalizedMag_dB,Nos,'moving');
% 
% % make sure column vector
% spectralTrend = spectralTrend(:);
% partOfNormalizedMag_dB = partOfNormalizedMag_dB(:);
% 
% % Subtract continuous part from partOfNormalizedMag_dB
% partOfNormalizedMag_dB = partOfNormalizedMag_dB-spectralTrend;
% 
% % Search the global maximum peak within range        
% [~, symratePeakIdx] = max(partOfNormalizedMag_dB);
% estimatedFineSymRate = partOfFreqAxies(symratePeakIdx);
% 
% % plot graph
% handles.data.x = partOfFreqAxies;
% handles.data.y = 10.^(partOfNormalizedMag_dB / 10);    
% handles.data.reqAxes = handles.axes_monitoring;
% handles.units = 'dB';
% plotGraph(handles);
% 
% hold(handles.axes_monitoring,'on');
% plot(partOfFreqAxies(symratePeakIdx), partOfNormalizedMag_dB(symratePeakIdx),'ro');
% % legend(handles.axes_monitoring, {'(2,1) cyclic moment','Symbol rate peak'});

%% save handles structure
handles.rx_sig = poweredSig;
guidata(hObject,handles);

%% update GUI states
set(hObject,'BackgroundColor',[222 235 250]/255);

function handles=pushbutton_coarseBandWidthEstimation_Callback(hObject, eventdata, handles)
%
% 주의 : 임계치보다 큰 스펙트럼에서 신호성분외의 것이 포함될때 오류가 난다.
%
updateLoggingMsg('While processing bandwidth estimation, please wait ...',handles);  
try
        
    %% execute spectrum callback    
    handles = pushbutton_viewSpectrum_Callback(hObject, eventdata, handles);
    %% Parse data from handles structure
    Pxx_dB              = handles.Pxx_dB;
    freqAxis            = handles.freqAxis;
    FFTlength           = handles.FFTlength;
    freqResolution      = handles.freqResolution;
    bandwidthConstant   = handles.bandwidthConstant;
    MINIMUM_HIST_COUNT  = handles.MINIMUM_HIST_COUNT;
    rx_sig = handles.rx_sig;

    % frequency smoothing
    smoothedPxx_dB      = smooth(Pxx_dB,1,'moving');  
        
    % Find min. and max. Duplicated code in histogram
    [yMax, yMaxIdx] = max(smoothedPxx_dB); yMin = min(smoothedPxx_dB);
    % Min-max normalization    

    MaxMinusMin = yMax-yMin;
%     binWidth = 3; % dB
    numOfBins = 35; % MaxMinusMin / binWidth; 
    binWidth = MaxMinusMin / numOfBins; % dB
    [nelements,xcenters_dB]=hist(smoothedPxx_dB,numOfBins);
%     figure; hist(smoothedPxx_dB,numOfBins);
   
    % Find noise floor
    [val, idx2] = max(nelements);
    noiseFloor = xcenters_dB(idx2);
    
    secondPeakCount = 0;
    offset = ceil(4/binWidth); % (offset * binWidth) dB
    threshold = (yMax*0.8 + yMin*0.2) ;
    if noiseFloor > threshold
        % In case of that power of signal is larger than power of noise 
        for i = 1:(idx2-4)
            % Check whether current peak is local maximum or not   
            
            % left index
            leftIdx = i-offset;  % (binWidth x offset) = ( ) dB
            if leftIdx<1, leftIdx= 1; end
            
            % right index
            rightIdx = i+offset; 
            if rightIdx>numOfBins, rightIdx = numOfBins; end
            
            window = leftIdx : rightIdx;   % from x-offset[dB] to x+offset[dB]
            isNotLocalMax = sum(nelements(window) > nelements(i));
            if nelements(i) > MINIMUM_HIST_COUNT
                if isNotLocalMax>0
                    % Skip current sample
                else
                    % Select current sample as second peak
                    if nelements(i) > secondPeakCount
                        secondPeakCount    = nelements(i);
                        secondPeakCenterdB = xcenters_dB(i);                
                    end                
                end           
            else
            end
        end
    else
        for i = numOfBins:-1:(idx2+4)
            % Check whether current peak is local maximum or not   
            leftIdx = i-offset;  if leftIdx<1, leftIdx= 1; end
            rightIdx = i+offset; if rightIdx>numOfBins, rightIdx = numOfBins; end
            window = leftIdx : rightIdx;
            if nelements(i) > MINIMUM_HIST_COUNT
                isNotLocalMax = sum(nelements(window) > nelements(i));

                if isNotLocalMax>0
                    % Skip current sample
                else
                    % Select current sample as second peak
                    if nelements(i) > secondPeakCount
                        secondPeakCount    = nelements(i);
                        secondPeakCenterdB = xcenters_dB(i);                
                    end     
                    break;
                end
            else
            end
        end
    
    end
    
    
    if secondPeakCount < 1
        % FM?
        [pxxMax, pxxMaxIdx] = max(xcenters_dB);
        signalFloorCount    = nelements(pxxMaxIdx);
        signalFloor         = pxxMax;
        
        bandwidthConstant = 0.3;
    else
        signalFloorCount    = secondPeakCount;
        signalFloor         = secondPeakCenterdB;       
    end
    
%     figure; hist(smoothedPxx_dB,numOfBins); 
%     xlabel('Magnitude of Spectrum [db]'); ylabel('Bin Count');
%     
%     hold on; plot(signalFloor,signalFloorCount,'ro');
%
    
%     numIdx = 1 : numOfBins;
%     th1 = 80/100*MaxMinusMin;% + yMin;
%     if noiseFloor > th1 % 신호 대역폭 similar to fs/2
%         signalFloor = noiseFloor;
%         levelB2NoiseNSig = signalFloor - offset;
%         
%         logicIdx = xcenters_dB < levelB2NoiseNSig;  
%         numIdx = numIdx(logicIdx);
%         [~, Idx]          = max(nelements(logicIdx));
%         locaMaxNumIdx     = numIdx(Idx);
%         
%         noiseFloor        = xcenters_dB(locaMaxNumIdx);
%         signalFloorCount  = nelements(locaMaxNumIdx);
%     else
%         logicIdx = xcenters_dB >= levelB2NoiseNSig;
%         numIdx = numIdx(logicIdx);
% 
%         [~, Idx]          = max(nelements(logicIdx));
%         locaMaxNumIdx     = numIdx(Idx);
%                 
%         signalFloor       = xcenters_dB(locaMaxNumIdx);
%         signalFloorCount  = nelements(locaMaxNumIdx);
%     end
    
    % find a minimum spectrum in the range of +-5 peak     
    const = 5;
    freamMin = min(smoothedPxx_dB(yMaxIdx-const: yMaxIdx+const)) ;
    lineSpectrumHeight = yMax-freamMin;
    
    % check whether sorrounding smaples from the peak are monotonic increase or decrease
    %  isMonotonic = 1;
    % C-style
%     for i= 1 : const 
%         tmpf =  (smoothedPxx_dB(yMaxIdx-const+i) - smoothedPxx_dB(yMaxIdx-const-1+i)); % left
%         tmpf2 = (smoothedPxx_dB(yMaxIdx+i) - smoothedPxx_dB(yMaxIdx+i-1));  % right
%         if tmpf < 0 || tmpf2 >0
%             isMonotonic = 0;
%         end
%     end 

    % check existence of another carrier components
    subLineSpecMagThreshold = yMax*0.7; % - lineSpectrumHeight/2; %dB by experimental 0.75*MaxMinusMin; % + yMin;
    isSingleLineSpectrum = 1; offset = 5;
    leftCon = sum(smoothedPxx_dB(1:yMaxIdx-offset) > subLineSpecMagThreshold);
    rightCon = sum(smoothedPxx_dB(yMaxIdx+offset:end) > subLineSpecMagThreshold);
    cond = leftCon + rightCon;
    if cond, isSingleLineSpectrum = 0; end
    
    th2 = 14; %dB , 140.3 * MaxMinusMin; % + yMin;
    if ( lineSpectrumHeight > th2 && isSingleLineSpectrum );  
        myMsg = 'Signal group having line-spectrum carrier component'; 
        existsCarrier = 1;        
 
        %% Get inst. amp.
        instAmp = abs(rx_sig); 
        
        instAmp=instAmp(10:end);        
        
        instAmp = smooth(instAmp,10,'moving');
        
        nBin = 35;
        [count, center]=hist(instAmp,nBin);
        [firstMaxCnt, firstMaxCntIdx] = max(count);
        
        isLeft2nd = 0;
        if (firstMaxCntIdx-1) > (nBin-firstMaxCntIdx)
           % Right normal dist.
           distWidth = nBin-firstMaxCntIdx;
%            idxVal = (firstMaxCntIdx-distWidth);

           idxVal = 2*firstMaxCntIdx-nBin;
           
           if idxVal < 1
               idxVal = 1;
           end
           
           range = 1 : idxVal;
           
           isLeft2nd = 1; % The bit indicates the location of 2nd dist. 
        else
           % Left normal dist.
           distWidth = firstMaxCntIdx-1;
%            idxVal = firstMaxCntIdx+distWidth;

           idxVal = 2*firstMaxCntIdx;
           
           if idxVal > nBin
               idxVal = nBin;
           end
           
           range = idxVal : nBin;
        end
        
        [secondMaxCnt,secondMaxCntIdx] = max(count(range));
        if isLeft2nd            
        else
            secondMaxCntIdx = secondMaxCntIdx + range(1);
        end
  
        peakDist = abs(firstMaxCntIdx-secondMaxCntIdx);
        if(0.5*firstMaxCnt < secondMaxCnt )% && (peakDist > distWidth)
            isAM = 0;
        else
            isAM = 1;
        end
        
%         if  signalFloorCount < 10  % AM            
%             isAM = 1;
%         else
%             offset = 14/100*MaxMinusMin; % + yMin;
%             idx=((signalFloor-offset) <= xcenters_dB )& (xcenters_dB <= (signalFloor+offset));
%             if sum( nelements(idx) > signalFloorCount) % AM
%                 sign = 1;
%             else % 2ASK
%                 sign = 0;
%             end
%             isAM = 0;
%         end

        if isAM % AM
            adFlag = 'Analog'; % analog or digital signal classification bits
            coarseModType = 'AM';
            noiseFloorOffset = 5/100*MaxMinusMin; % + yMin;
            bwThreshold_dB = noiseFloor + noiseFloorOffset; % no reason                        
            
        else % 2ASK
            adFlag = 'Digital'; % analog or digital signal classification bits
            coarseModType = '2ASK';
            bwThreshold_dB = (signalFloor - noiseFloor) * ...
                bandwidthConstant + noiseFloor;
        end
    else         
       % others           
        myMsg = 'Signal group not having line-spectrum carrier component';
        existsCarrier = 0;
        
        coarseModType = 'unKnown';
        adFlag = 'Digital'; % analog or digital signal classification bits
        bwThreshold_dB = (signalFloor - noiseFloor) * ... 
            bandwidthConstant + noiseFloor;
    end

    %% convert spectrum to linear scale 
    smoothedPxx = 10.^(smoothedPxx_dB / 10);
    bwThreshold   = 10.^((bwThreshold_dB) / 10);

    num = 1: FFTlength;
    idx =  smoothedPxx >= bwThreshold;
    freqIdx = num(idx);
    binMinSize = 1400 / freqResolution;
    [~,yMaxIdx] = max(smoothedPxx(idx));
    % 대역폭내 주파수 빈 최대 사이 간격 500Hz/33Hz  = 151      
    % 원치않은  스펙트럼 성분이 다수라면 중간값부터 밖으로 나가면서 모두 제거한다.
    
    % toward left edge
    partIdx = freqIdx(1:yMaxIdx); 
    subIdx = partIdx(2:end) - partIdx(1:end-1);
    num = 1 : length(partIdx);
    cond1 = num(subIdx > binMinSize); 
    if ~isempty(cond1), axisLowerIdx = partIdx(max(cond1)+1);
    else axisLowerIdx = partIdx(1);
    end
    
    % toward right edge
    partIdx = freqIdx(yMaxIdx:end); 
    subIdx = partIdx(2:end) - partIdx(1:end-1);
    num = 1 : length(partIdx);
    cond1 = num(subIdx > binMinSize); 
    if ~isempty(cond1), axisUpperIdx =partIdx(min(cond1));
    else axisUpperIdx = partIdx(end);   
    end
  
    updateLoggingMsg(['Modulation type is ', coarseModType],handles) 
    updateLoggingMsg(myMsg,handles)
    
    BWSampleIdx = axisLowerIdx : axisUpperIdx;
    
    if strcmp(coarseModType,'AM')
        coarseBandWidth = 15000;
    else
        coarseBandWidth = freqAxis(axisUpperIdx) - freqAxis(axisLowerIdx);
    end
    
    %% SNR estimation
%     minB = BWSampleIdx(1) - 3;  maxB = BWSampleIdx(end) + 3;
%     BWSampleIdx = minB : maxB;
    noiseSampIedx                  = 1 : FFTlength;
    noiseSampIedx(BWSampleIdx) = [];

    signalPower = mean( db2pow(Pxx_dB(BWSampleIdx)) );
    noisePower = mean( db2pow(Pxx_dB(noiseSampIedx)) );

    estimatedSNR_dB    = (10.0)*log10(signalPower / noisePower);
    
if ~handles.opMode
    %% update GUI states
    set(handles.pushbutton_coarseSNRestimation,'Enable','on');
    set(handles.pushbutton_addNoise,'Enable','on');
    set(handles.pushbutton_applyFading,'Enable','on');
   
    set(hObject,'BackgroundColor',[222 235 250]/255);
    str1 = ['Estimated bandwidth : ',num2str(coarseBandWidth), ' Hz'];
    updateLoggingMsg(str1,handles);
     
    %% plot graph
    cla(handles.axes_monitoring,'reset');
    legend(handles.axes_monitoring,'off');

    thLine_dB = ones(1,length(freqAxis))*bwThreshold_dB;

    plot(freqAxis,smoothedPxx_dB,'parent',handles.axes_monitoring);
    axis(handles.axes_monitoring,[min(freqAxis), max(freqAxis), min(smoothedPxx_dB), max(smoothedPxx_dB)]);
    ylabel('Normalized magnitude','parent',handles.axes_monitoring);
    xlabel('Frequency (Hz)','parent',handles.axes_monitoring);  
    
    hold(handles.axes_monitoring, 'on'); grid(handles.axes_monitoring, 'on');
    % horizontal threshold line
    plot(freqAxis,thLine_dB,'-r','parent',handles.axes_monitoring);
    % vertical threshold line
    plot(ones(length(smoothedPxx_dB),1)*freqAxis(axisUpperIdx) ...
        ,smoothedPxx_dB,'-r','parent',handles.axes_monitoring);
    plot(ones(length(smoothedPxx_dB),1)*freqAxis(axisLowerIdx) ...
        ,smoothedPxx_dB,'-r','parent',handles.axes_monitoring);
    % samples in BW
    plot(freqAxis(BWSampleIdx),smoothedPxx_dB(BWSampleIdx),'y.','parent',handles.axes_monitoring);
    plot(freqAxis(axisUpperIdx),thLine_dB(axisUpperIdx),'go','parent',handles.axes_monitoring);
    plot(freqAxis(axisLowerIdx),thLine_dB(axisLowerIdx),'k*','parent',handles.axes_monitoring);
    % legend(handles.axes_monitoring, {'Smoothed PSD','Horizontal threshold', ...
    %     'Vertical threshold','Samples larger than threshold'});
    updateLoggingMsg(['Estimated SNR : ', num2str(estimatedSNR_dB),' dB'],handles);
    updateLoggingMsg('Complete bandwidth estimation',handles); 
end

    %% save data
    handles.coarseBandWidth      = coarseBandWidth ;
    handles.BWSampleIdx          = BWSampleIdx;
    handles.smoothedPxx_dB       = smoothedPxx_dB;
    handles.freqAxis             = freqAxis;
    handles.BWthreshold          = bwThreshold_dB;
    handles.adFlag               = adFlag;
    handles.existsCarrier        = existsCarrier;
    handles.coarseModType        = coarseModType;
    handles.signalFloor          = signalFloor;
    handles.noiseFloor           = noiseFloor;
    handles.estimatedSNR         = estimatedSNR_dB;  
    handles.signalPower          = signalPower;
    handles.noisePower           = noisePower;
    
    guidata(hObject,handles);  %% method  1 - x-dB point bandwidth measurement
    % [maxPxxdB, idxMax] = max();
    % updateLoggingMsg('Choose the proper x-dB point',handles);  
    % [x,currentYpointOnAxes] = ginput(1);
    % XdBPoint = maxPxxdB - currentYpointOnAxes;
    % threshold_dB = currentYpointOnAxes;
    % set(handles.edit_xdB,'String',num2str(XdBPoint));
    % threshold          = 10^(threshold_dB/10);
    % 
    % num = 1: FFTlength;
    % tmp = num( smoothedPxx_dB >= threshold_dB );
    % 
    % axisLowwerIdx = tmp(1);
    % axisUpperIdx = tmp(end);
    % 
    % BWSampleIdx = axisLowwerIdx : axisUpperIdx;
    % coarseBandWidth = freqAxis(axisUpperIdx) - freqAxis(axisLowwerIdx);

    %% method 1.5
    % segIndex = 1 : FFTlength;
    % idxHigherThanThreshold = segIndex(BWSampleIdx);
    % % divide each spectrum region 
    % diffIdx = diff(idxHigherThanThreshold);
    % idxOfIdx = diffIdx > 20; % important constant
    % axisLowwerIdx =[ idxHigherThanThreshold(idxOfIdx) idxHigherThanThreshold(end)];
    % IdxOfIdx2 = logical([ 0 idxOfIdx(1:end-1)]);         
    % axisUpperIdx =[ idxHigherThanThreshold(1) idxHigherThanThreshold(IdxOfIdx2)]; 
    % 
    % % exception handling
    % rightlength = length(axisUpperIdx); leftlength = length(axisLowwerIdx);
    % if rightlength ~= leftlength
    %     minlength = min(leftlength,rightlength);
    %     axisUpperIdx = axisUpperIdx(1:minlength);
    %     axisLowwerIdx = axisLowwerIdx(1:minlength);
    % end
    % 
    % coarseBandWidth = sum( abs( freqAxis(axisUpperIdx) - freqAxis(axisLowwerIdx) ));

       %% method 2
        % get lowwer bound
%         [~, idxMax] = max(smoothedPxx_dB);
%         for jj = 1 : idxMax
%             if ( smoothedPxx_dB( idxMax -jj ) < threshold_dB )
%                 axisLowerIdx = idxMax -jj ;
%                 break;
%             end            
%         end        
%         % get upper bound
%         for kk = (idxMax+1) : FFTlength
%             if (smoothedPxx_dB(kk) < threshold_dB)
%                 axisUpperIdx = kk;
%                 break;
%             end
%         end   
%         axisLowerIdx = tmp(1);
%         axisUpperIdx = tmp(end);

    %% Method3 : integration-based bandwidth measurement
    % sorted_smoothedPxx = sort(smoothedPxx,'descend');
    % integPxx = cumsum(sorted_smoothedPxx);
    % percentage = integPxx / integPxx(end) * 100;
    % 
    % % for i = 1 : length(percentage)
    % %     if percentage(i) >50
    % %         beginIdx = i; break;
    % %     end
    % % end
    % 
    % beginIdx = 1;
    % % line equation
    % x1=percentage(1);            x2=percentage(end);
    % y1=sorted_smoothedPxx(1);    y2= sorted_smoothedPxx(end);
    % slope = (y2-y1) / (x2-x1);
    % ycoef = 1;
    % xcoef = -slope;
    % const = -y1 + slope*x1;
    % 
    % % distance between two points.            
    % dist = (xcoef*percentage(beginIdx:end) + ycoef*sorted_smoothedPxx(beginIdx:end) + const )...
    %     / ( sqrt( xcoef^2 +ycoef^2));
    % 
    % % find a point 
    % [minVal, minIdx] = min(dist);
    % 
    % threshold = sorted_smoothedPxx(minIdx -5);
    % threshold_percentage = percentage(minIdx);
    % 
    % BWSampleIdx = smoothedPxx >= threshold;
    % 
    % segIndex = 1 : FFTlength;
    % idxHigherThanThreshold = segIndex(BWSampleIdx);
    % % divide each spectrum region 
    % diffIdx = diff(idxHigherThanThreshold);
    % idxOfIdx = diffIdx > 20; % important constant
    % axisLowwerIdx =[ idxHigherThanThreshold(idxOfIdx) idxHigherThanThreshold(end)];
    % IdxOfIdx2 = logical([ 0 idxOfIdx(1:end-1)]);         
    % axisUpperIdx =[ idxHigherThanThreshold(1) idxHigherThanThreshold(IdxOfIdx2)]; 
    % 
    % % exception handling
    % rightlength = length(axisUpperIdx); leftlength = length(axisLowwerIdx);
    % if rightlength ~= leftlength
    %     minlength = min(leftlength,rightlength);
    %     axisUpperIdx = axisUpperIdx(1:minlength);
    %     axisLowwerIdx = axisLowwerIdx(1:minlength);
    % end
    %             
    % coarseBandWidth = sum( abs( freqAxis(axisUpperIdx) - freqAxis(axisLowwerIdx) ));

    %% method 4
    % numOfBins = 50;
    % [nelements,xcenters]=hist(smoothedPxx_dB,numOfBins);
    % 
    % yMax = max(smoothedPxx_dB); yMin = min(smoothedPxx_dB); 
    % c1 = 0.4;    
    % threshold = yMin + c1 * (yMax-yMin);

    % find local maximum
    % logicIdx = xcenters > threshold;
    % numIdx = 1 : numOfBins; numIdx = numIdx(logicIdx);
    % [val, Idx] = max(nelements(logicIdx));
    % locaMaxNumIdx = numIdx(Idx);
    % nelements(locaMaxNumIdx)
    % localMaximumSpectrum = xcenters(locaMaxNumIdx);
    % delta =  nelements(locaMaxNumIdx) - median(nelements(logicIdx));  
    % 
    % % rectangular pulse
    % numOfSamplesOfminBW = minimumSymbolRate / freqResolution;
    % const = 0.8;
    % if  delta < floor(const*numOfSamplesOfminBW)
    %     X = sprintf('Analog modulation signal, %d',delta);        
    % else
    %     X = sprintf('Digital modulation signal, %d',delta);        
    % end  
    % disp(X);
    % threshold = localMaximumSpectrum;
   
catch err
    updateLoggingMsg('Error while coarse bandwidth estimation section ',handles);
    rethrow(err);    
end

function handles=pushbutton_coarseSNRestimation_Callback(hObject, eventdata, handles)
updateLoggingMsg('While processing SNR estimation, please wait ...', handles);
try
    signalSampleIdx                = handles.BWSampleIdx;
%     minB = signalSampleIdx(1) - 3;
%     maxB = signalSampleIdx(end) + 3;
%     signalSampleIdx = minB : maxB;
    Pxx_dB                         = handles.Pxx_dB;
    FFTlength                      = handles.FFTlength;
    noiseSampIedx                  = 1 : FFTlength;
    noiseSampIedx(signalSampleIdx) = [];
    freqAxis                       = handles.freqAxis;
%     featureSNR_dB                  = handles.featureSNR_dB;  
    %%% in case of intergration method
    % sigPower_dB = max( smoothedPxx_dB(etimatedSampleIdx) );
    % noisePower_dB = mean( smoothedPxx_dB(~etimatedSampleIdx));
    sigPower = mean( db2pow(Pxx_dB(signalSampleIdx)) );
    noisePower = mean( db2pow(Pxx_dB(noiseSampIedx)) );
    
%     sigPower_dB = mean( Pxx_dB(signalSampleIdx) );
    %noisePower_dB = mean( Pxx_dB(noiseSampIedx));

    estimatedSNR_dB    = (10.0)*log10(sigPower / noisePower);
%     estimatedSNR_dB = sigPower_dB - noisePower_dB;

    %% convert estimatedSNR to index
%     dd = ( featureSNR_dB(1:end-1) + featureSNR_dB(2:end) ) /2;
%     dd(end+1) = inf;
%     for i= 1 : length(dd)
%         if estimatedSNR_dB < dd(i)
%              estimatedSNRIdx = i;   break;
%   _      end
%     end
    
if ~handles.opMode    % execute only manual mode
    %% update GUI states
    set(handles.pushbutton_coarseSNRestimation,'BackgroundColor',[222 235 250]/255);
    
    updateLoggingMsg(['Estimated SNR : ', num2str(estimatedSNR_dB),' dB'],handles);
    %% plot 
    cla(handles.axes_monitoring,'reset');
    legend(handles.axes_monitoring,'off');
    plot(freqAxis,Pxx_dB,'parent',handles.axes_monitoring);
    hold(handles.axes_monitoring,'on'); grid(handles.axes_monitoring,'on');
    plot(freqAxis(signalSampleIdx), Pxx_dB(signalSampleIdx),'g.','parent',handles.axes_monitoring);
    plot(freqAxis(noiseSampIedx), Pxx_dB(noiseSampIedx),'rs','parent',handles.axes_monitoring);
    xlim(handles.axes_monitoring, [min(freqAxis), max(freqAxis)]);  
    ylim(handles.axes_monitoring,[min(Pxx_dB), max(Pxx_dB)]);
    legend(handles.axes_monitoring,{'Spectrum','Signal samples','Noise samples'});
    ylabel(handles.axes_monitoring,'Magnitude(dB)');
    hold(handles.axes_monitoring,'off');
    
    set(handles.pushbutton_coarsefreqOffsetEstimation,'Enable','on');
end
    %% update GUI handles
    guidata(hObject,handles);

%      handles.sigPower_dB     = floor(sigPower_dB);
    handles.estimatedSNR    = estimatedSNR_dB;
%     handles.estimatedSNRIdx = estimatedSNRIdx;
    guidata(hObject,handles);
catch err
    updateLoggingMsg('Error while processing SNR estimation'.',handles);
    rethrow(err);
end

function handles=pushbutton_addNoise_Callback(hObject, eventdata, handles)
updateLoggingMsg('While processing SNR adjustment, please wait ...',handles);

try
    SNR_dB       = handles.SNR_dB;
    rx_sig       = handles.basic_rx_sig;
%     sam_freq     = handles.sam_freq;
%     signalFloor  = handles.signalFloor;
%     noiseFloor   = handles.noiseFloor;
    signalPower  = handles.signalPower;
%     noisePower   = hnadles.noisePower;
    
    viewEndIdx   = length(rx_sig); 
    viewStartIdx = 1 ;
    displayRange = viewStartIdx:viewEndIdx;    
    % make average power to 1
    signalPower = 10*log10(signalPower);
    noisy_rx_sig = awgn(rx_sig(displayRange),SNR_dB,signalPower);
    
    %% save handles structures
    handles.rx_sig       = noisy_rx_sig;
    handles.viewStartIdx = viewStartIdx;
    handles.viewEndIdx   = viewEndIdx;
    guidata(hObject,handles);
    if ~handles.opMode

        %% update GUI states
        cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off');
       % execute spectrum callback    
        handles = pushbutton_viewSpectrum_Callback(hObject, eventdata, handles);
        
        % set(handles.pushbutton_preLPF,'Enable','on');
        set(handles.edit_viewStartIdx,'String',num2str(viewStartIdx));
        set(handles.edit_viewEndIdx,'String',num2str(viewEndIdx));
        set(handles.edit_inputNumSamples,'String',num2str(viewEndIdx));
     
        set(handles.pushbutton_coarsefreqOffsetEstimation,'Enable','on');
        set(handles.pushbutton_coarseSNRestimation,'Enable','on');
        
        set(hObject,'BackgroundColor',[222 235 250]/255);
        updateLoggingMsg('Complete SNR adjustment',handles);
    end
catch err
    updateLoggingMsg('Error in SNR adjustment',handles);
    rethrow(err);    
end

function slider_snr_CreateFcn(hObject, ~, handles)
set(hObject,'Value',40);    
set(hObject,'Min',-20);  
set(hObject,'Max',40);
minorStep = 0.016;
majorStep = 0.1;
set(hObject,'SliderStep',[minorStep, majorStep]);

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function handles=slider_snr_Callback(hObject, eventdata, handles)
snr_dB       = floor( get(hObject,'Value') );
%% update GUI states
set(handles.edit_dispSNR,'String',num2str(snr_dB));
%% save data
handles.SNR_dB = snr_dB;
guidata(hObject, handles);

function edit_fadingCarrierFreq_Callback(hObject, eventdata, handles)
fo = 1e6*str2double( get(hObject,'String') ); % in MHz
v  = handles.averageSpeedOfObject;
c  = 10.8e8;                           % The speed of light  [km/h]   
handles.maxDopplerShift  = floor(( fo * v) / c);
handles.fadingCarierFreq = fo;

set(handles.edit_maxDopplerShift,'String',num2str(handles.maxDopplerShift));
guidata(hObject,handles);

function edit_fadingCarrierFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_fadingObjSpeed_Callback(hObject, eventdata, handles)
v = str2double( get(hObject,'String') );
c = 10.8e8;                        % The speed of light  [km/h]    
handles.maxDopplerShift = (handles.fadingCarierFreq * v) / c;
handles.averageSpeedOfObject  = v;
set(handles.edit_maxDopplerShift,'String',num2str(handles.maxDopplerShift));
guidata(hObject,handles);

function edit_fadingObjSpeed_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_maxDopplerShift_Callback(hObject, eventdata, handles)
handles.maxDopplerShift= eval( get(hObject,'String') );
guidata(hObject,handles);

function edit_maxDopplerShift_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_fadingType_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String')); % returns popupmenu contents as cell array
fadingType = (contents{get(hObject,'Value')}); % returns selected item from popupmenu
handles.fadingType = fadingType;
guidata(hObject,handles);

function popupmenu_fadingType_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton_applyFading_Callback(hObject, ~, handles)
updateLoggingMsg('While processing fading channel, please wait ...',handles);
try
    ts           = 1/ handles.sam_freq;
    fd           = handles.maxDopplerShift;     % v*fo/c;    % Maximum Doppler shift [Hz]
    rx_sig       = handles.rx_sig;
    viewEndIdx   = handles.viewEndIdx;
    viewStartIdx = handles.viewStartIdx;
    %% the maximum Doppler shift      corresponding to vehicle speeds   at 900 MHz carrier frequency                
    % 27Hz                                  20 mph =  8.94 m/s =  32.186 km/h
    % 67Hz                                  50 mph = 22.35 m/s =  80.46 km/h
    % 107Hz                                 80 mph = 36.76 m/s = 128.74 km/h

    displayRange = viewStartIdx:viewEndIdx;   
    rx_sig = rx_sig(displayRange);
     % average path gains
%     c1.DopplerSpectrum = doppler.flat;
    c1.NormalizePathGains = true; 
    switch handles.fadingType 
        case 'Rayleigh'
           c1= rayleighchan(ts,fd);  
           c1.PathDelays =[0 ]; % Change th number of delays
           % Matlab automatically changed the size of c1.AvgPathGaindB,
           % c1.PathGains, and c1.ChannelFilterDelay
            fadedSig = filter(c1, rx_sig);
        case 'Rician'
           c1= ricianchan(ts,fd,0);
           
            fadedSig = filter(c1, rx_sig);
        case 'None'
           fadedSig = rx_sig;
        otherwise
    end
       
   
    % normalize signal power to 1
%     fadedSig = sqrt(length(fadedSig)) * fadedSig / norm(fadedSig,2);       
    %% save handles structures
    handles.rx_sig          = fadedSig;
    guidata(hObject,handles);
    
     if ~handles.opMode
        %% update GUI states
        cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off');
        set(handles.edit_inputNumSamples,'String',num2str(viewEndIdx));
        set(hObject,'BackgroundColor',[222 235 250]/255);
        updateLoggingMsg('Complete fading effect',handles);
     end
    
catch err
    updateLoggingMsg('Error in Fading channel',handles);
    rethrow(err);    
end

function handles=pushbutton_coarsefreqOffsetEstimation_Callback(hObject, eventdata, handles)
updateLoggingMsg('While processing frequency offset estimation , please wait ...',handles);
try
    rx_sig                  = handles.rx_sig;
    Pxx_dB                  = handles.Pxx_dB;
    freqAxis                = handles.freqAxis ;
    BWSampleIdx             = handles.BWSampleIdx;
    FFTlength               = handles.FFTlength;
    existsCarrier           = handles.existsCarrier;
    rx_sig = rx_sig(:);     

    if existsCarrier 
        % Based on maximum point
        [~, offsetIdx] = max(Pxx_dB); 
        carrierFreqOffset = freqAxis(offsetIdx);        
    else
       % Based on middle point of esimtaed BW                  
        offsetIdx = floor(median(BWSampleIdx));
        carrierFreqOffset = freqAxis(offsetIdx);        
    end
        
    bwShiftIdx = offsetIdx - (FFTlength/2 +1);
    BWSampleIdx = BWSampleIdx - bwShiftIdx;
    
    % Shift spectrum toward zero frequency
    ttIdx = (0 : length(rx_sig)-1 );
    rx_sig = rx_sig.*exp( -1i*2*pi*bwShiftIdx/FFTlength*ttIdx ).';

    
    %% save handle structure
    handles.rx_sig           = rx_sig;
    handles.coarseFreqOffset = carrierFreqOffset;    
    handles.BWSampleIdx      = BWSampleIdx;
    guidata(hObject,handles);
    
    % Get new spectrum
    handles = pushbutton_viewSpectrum_Callback(hObject, eventdata, handles);      
    
%     figure; plot(handles.freqAxis,handles.Pxx_dB);
    if ~handles.opMode       
        %% update PSD plot
%         cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off'); hold(handles.axes_monitoring,'on');

        plot(freqAxis,Pxx_dB,'r--','parent',handles.axes_monitoring); % previous spectrum
        axis(handles.axes_monitoring,[min(freqAxis), max(freqAxis), min(Pxx_dB),max(Pxx_dB)]);
        hold(handles.axes_monitoring,'on'); grid(handles.axes_monitoring,'on');

        %% Plot lines
        plot(ones(length(Pxx_dB),1)*freqAxis(offsetIdx) ...
            ,Pxx_dB,'-r','parent',handles.axes_monitoring);
        plot(ones(length(Pxx_dB),1)*freqAxis(floor(end/2)+1) ...
            ,Pxx_dB,'-.g','parent',handles.axes_monitoring);
        
%         legend(handles.axes_monitoring,{'Spectrum','Offset frequency','Center frequency'});      
        ylabel(handles.axes_monitoring,'Magnitude(dB)');
        xlabel(handles.axes_monitoring,'Frequency (Hz)');
        hold(handles.axes_monitoring,'off');

        %% update GUI states
        set(handles.pushbutton_lowPassFilterling,'Enable','on');
        set(hObject,'BackgroundColor',[222 235 250]/255);
        updateLoggingMsg(['Estimated carrier frequency offset : ', ...
            num2str(carrierFreqOffset),' Hz'],handles);
    end
catch err
    updateLoggingMsg('Error while processing frequency offset estimation ',handles);
    rethrow(err);    
end

function handles = pushbutton_lowPassFilterling_Callback(hObject, ~, handles)
updateLoggingMsg('While processing lowpass filtering, please wait ...',handles);
try   
    % 반드시 신호의 주파수 대역만 추가 되어야 한다.
    sam_freq                = handles.sam_freq;
    BWSampleIdx             = handles.BWSampleIdx;
    freqAxis                = handles.freqAxis;     
    adFlag                  = handles.adFlag;    
    rx_sig                  = handles.rx_sig;
    cutOffExpansionRatio    = handles.cutOffExpansionRatio;    
    
    %% baseband filters
    bwFreqAxis = freqAxis(BWSampleIdx);
%     lowerBound = bwFreqAxis(1);
    upperBound = bwFreqAxis(end);
   
    switch adFlag
        case 'NoSignal'
            filteredSig = rx_sig;
        case {'Analog','Digital'}    
            filterOrder = 31;  
            Wn = abs ( cutOffExpansionRatio*upperBound / (sam_freq/2 ) ) ;
          
            lowPassFilterCoefficients = fir1(filterOrder,Wn);
            lowPassFilterCoefficients = lowPassFilterCoefficients(:);

            % magnitudeSpectrumOfFilter=fftshift(fft(lowPassFilterCoefficients));
            % lowPassFilterCoefficients = lowPassFilterCoefficientswindow
            %% filtering method using filter function in MATLAB
            filteredSig = filter(lowPassFilterCoefficients,1,rx_sig);
            
%             re = filter(lowPassFilterCoefficients,1,real(rx_sig));
%             im = filter(lowPassFilterCoefficients,1,imag(rx_sig));
%             fSig = re+1i*im;
        otherwise
    end


    %% filtering method using convolution
    % numberOfSegments = (length(rx_sig))./(FFTlength);
    % numberOfSegments = fix(numberOfSegments);
    % 
    % xStart = 1:FFTlength:numberOfSegments*FFTlength;
    % xEnd   = xStart+FFTlength-1;
    % 
    % filteredSig = []; % zeros(1 : numberOfSegments * FFTlength);
    % 
    % for i = 1:numberOfSegments
    % %     win = hanning(length(FFTlength));
    %     xw = rx_sig(xStart(i):xEnd(i));
    %     
    %     % method 1
    % %     filteredSigPart = conv(xw,lowPassFilterCoefficients);  
    % %     filteredSig = [filteredSig; filteredSigPart];
    % 
    %     % method 2
    %     filteredSigPart = conv(xw,lowPassFilterCoefficients,'same');  
    %     filteredSig(xStart(i):xEnd(i) ) = [filteredSigPart];      % Auto spectrum.
    % end
    % figure; scatterplot(filteredSig);
    % figure; plot(real(filteredSig));

    if ~handles.opMode
        %% update PSD plot
        cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off');

        %% update GUI states
        set(handles.pushbutton_coarseSymRateEstimation,'Enable','on');
        
        set(hObject,'BackgroundColor',[222 235 250]/255);
        updateLoggingMsg('Complete lowpass filtering',handles);
    end
    %% update handles structure
    handles.rx_sig = filteredSig;
    guidata(hObject,handles);
catch err
    updateLoggingMsg('Error in filtering',handles);
    rethrow(err);    
end

function handles = pushbutton_coarseSymRateEstimation_Callback(hObject, ~, handles)
    updateLoggingMsg('While processing coarse symbol rate estimation , please wait ...',handles);
try
    % Parse the data
    rx_sig               = handles.rx_sig;
    sam_freq             = handles.sam_freq;
    minimumSymbolRate    = handles.minimumSymbolRate;    
    maximumSymbolRate    = handles.maximumSymbolRate;   
    viewEndIdx           = handles.viewEndIdx;
    viewStartIdx         = handles.viewStartIdx;
    NFFT                 = handles.FFTlength;
    estBW                = handles.coarseBandWidth;
    adFlag               = handles.adFlag;
    coarseModType        = handles.coarseModType;
    
    sigmaAThreshold     = handles.sigmaAThreshold;
    lowerConstant       = handles.lowerConstantBW;  
    upperConstant       = handles.upperConstantBW; 
    maxixumBandwidthForLinearModulation = handles.maxixumBandwidthForLinearModulation;
    isDehoppingSig      = handles.isDehoppingSig;
    
    smallSamplesPerFrame = 4000;
    if isDehoppingSig 
        % 역도약 신호는 FM으로 가정
        coarseModType = 'FM';                     
        coarseSamplesPerSymbol  = 0;
        coarseNumberOfsymbols   = 0;
        coarseSymbolRate        = 0;  
        sigma_a                 = 0;
        linearFlag             = 'Nonlinear';    % linear digital signal classification bits       
        adFlag = 'Analog';        
        updateLoggingMsg('Modulation type is FM',handles);       
        
    else
        % Restrict range of signal
        rx_sig              = rx_sig(viewStartIdx:viewEndIdx);
        if length(rx_sig) < smallSamplesPerFrame
            smallSamplesPerFrame = length(rx_sig); numSegment=1;
        else
            numSegment = floor( length(rx_sig) / smallSamplesPerFrame ); 
        end   

        if numSegment > 10, numSegment = 10; end % 15.06.09
        
       %% Classification between amplitude and angle modulation based on instantaneous amplitude        
        sigma_a = zeros(1,numSegment);      % Preallocate memory
        segIndex = 1 : smallSamplesPerFrame;
        rmsDifference = zeros(1,numSegment);
        for frameCount=1:numSegment
            % Divide signal with size of samplesPerFrame
            partialSig = rx_sig(segIndex);
            if frameCount == 1
                tmp = partialSig(end);
                partialSig = partialSig(50:end);
                partialSig(end+1:end+49) = tmp;
            end

            
    %         figure; plot(abs(partialSig)); pause(0.5); close(gcf);
            segIndex = segIndex+smallSamplesPerFrame;          
            % instantaneous amplitude
            a_amplitude = abs(partialSig); 
            
            a_amplitude= smooth(a_amplitude,3,'moving');
            
            m_a=mean(a_amplitude);                        % Eq. (2.3) in [1]
            a_n=a_amplitude/m_a;                          % Eq. (2.2) in [1]
            a_cn=a_n-1;                                   % Eq. (2.2) in [1]              
    %         % 6. sigma_a
    %         % Standard deviation of the normalized-centered instantaneous in the non-weak intervals of a signal segment
            sigma_a(frameCount)=sqrt(mean(a_cn.^2)-(mean(a_cn))^2);   % Eq. (4.1) at PP.118 in [1]
%             figure; plot(a_amplitude)
            y = a_amplitude ;
            x = 1 : smallSamplesPerFrame;
            v =[x.' , ones(smallSamplesPerFrame,1)];
            % increasing order
            p = v\y;    % refer to iLinearFit.m
            % Same as 
            % [Q,R] = qr(v,0);
            % p2 = full(R\(Q'*y));    % Same as p = D*A\(D*y);
            haty = v*p;
            rmsDifference(frameCount) =  sqrt(mean((y - haty).^2) );  
        end
        sigma_a = sigma_a(:);
        sigma_a = dataConditioning(sigma_a);  %figure; boxplot(sigma_a);
    %       rmsDifference =dataConditioning(rmsDifference);figure; boxplot(sigma_a);
       %% Classify between linear and non-linear modulation signal
    %    sigma_a = 0.001;
        if  median(sigma_a) > sigmaAThreshold
            linearFlag = 'Linear'; 
            coarseModType = 'unKnown';
            updateLoggingMsg('Linear modulation signals (PSK, QAM, PAM(ASK))',handles); 

             if estBW > maxixumBandwidthForLinearModulation
                 % How to handles this case?
                  linearFlag = 'nonLinear'; 
                  updateLoggingMsg('Linear modulation signals exceeding maximum bandwidth --> nonlinear',handles);
             end
        else
            linearFlag = 'nonLinear'; 
            updateLoggingMsg('Non-linear modulation signals',handles);
        end
        
       %% Symbol rate estimation 
        if length(rx_sig) < NFFT
            NFFT = 2.^nextpow2(length(rx_sig)); numSegment=1;
            segIndex = 1 : length(rx_sig);
        else
            numSegment = floor( length(rx_sig) / NFFT ); 
            segIndex = 1 : NFFT;
            if numSegment > 19, numSegment = 19; end
        end     
        
        % Preallocate memory    
        averagedCM21 = 0; % zeros(samplesPerFrame,1); 
        freqAxis = (0:NFFT/2-1) /NFFT * sam_freq ; 
        
        for frameCount=1:numSegment
            partialSig = rx_sig(segIndex);
            segIndex = segIndex+NFFT; 
            if  strcmp(linearFlag,'nonLinear');                % divide signal with size of samplesPerFrame
                % Obtain instantaneous frequency
                phi_phase = angle(partialSig);  % between +- pi
                
%                 figure; plot(phi_phase,'--','Marker','.');
                
                iFreq = diff( unwrap( phi_phase ) );   
                % noise reduction
                iFreq = smooth(iFreq,25,'moving');  
                if frameCount ==1
                % value of front some samples is incorrect. remove them by shifting samples. 
                % tail samples are replaced with a last sample  
                % 50 -> 70 at 15.06.16
                    iFreq(1:end-70)=  iFreq(71:end);
                    iFreq(end-71:end)= iFreq(end); 
                end
                 partialSig = iFreq;
%                  figure; plot( iFreq);
            end            
            %partialSig = partialSig - mean(partialSig);    
           
           %% Oerder and H.Meyr algorithm = Squaring method in MATLAB
            nonLinearTransformedX = partialSig.* conj(partialSig);

            Xk = 1/NFFT * fft(nonLinearTransformedX, NFFT); % fourier coefficient
            CM21 =  abs(Xk);  %figure; plot(10*log10(CM21));
            averagedCM21 = averagedCM21 + CM21;
        end
        averagedCM21 = averagedCM21 ./ numSegment;
        averagedCM21_dB = 10*log10(averagedCM21);
    
       %% Search symbol rate in range from min to max symbol rate
        % C style coding  , fs >= 4*BW should be satisfied !!!, 
        count=0; partialCM21_dB=[]; partialFreqAxis=[];
        if  strcmp(linearFlag,'nonLinear');                        
            % Move spectrum near to zero and search peaks in the range from
            % minimum symbol rate and maximum symbol rate
            for i = 1 : NFFT/2
               if minimumSymbolRate <= freqAxis(i) && ...
                       freqAxis(i) <= maximumSymbolRate
                   %if count == 0, offset= i-1;   end
                   count = count + 1;
                   partialCM21_dB(count) = averagedCM21_dB(i);
                   partialFreqAxis(count) = freqAxis(i);       
               end    
            end
        else
            % Move spectrum near to zero            
            for i = 1 : NFFT/2
               if  lowerConstant*estBW <= freqAxis(i) && ...
                       freqAxis(i) <= upperConstant*estBW 
                   %if count == 0, offset= i-1;   end
                   count = count + 1;
                   partialCM21_dB(count) = averagedCM21_dB(i);
                   partialFreqAxis(count) = freqAxis(i);       
               end    
            end
        end
%         figure; plot(freqAxis,averagedCM21_dB(1:end/2));
%         xlabel('Frequency (Hz)');
%         ylabel('Magnitude (dB)');
        % use a moving average filter with a Nos span to smooth all the data at once
        Nos = 13; %round( count * 0.5); %filter span        
        spectralTrend = smooth(partialCM21_dB,Nos,'moving');
 
        % make sure column vector
        spectralTrend = spectralTrend(:);
        partialCM21_dB = partialCM21_dB(:);

        % Subtract continuous part from partOfNormalizedMag_dB
        detrendedPartialCM21_dB = (partialCM21_dB)-spectralTrend; % Detrend data
%         figure; plot( partialCM21_dB);
%         figure; plot(abs(diff(partialCM21_dB))); title('Diff. spectrum');
        % Search the global maximum peak within range        
        [symratePeak, symratePeakIdx] = max(detrendedPartialCM21_dB);
        if  strcmp(linearFlag,'nonLinear')
            % Check harmonic symbol rate line spectrum for rectangular
            % pulse 
            isRectPulse= 0 ;          
            MINPEAKDIST = floor(minimumSymbolRate / handles.freqResolution);           
            AMP_THRESHOLD = symratePeak / 2;              
            [NumPeaks, sortedLocs]= ...
                 peakDetection(detrendedPartialCM21_dB,MINPEAKDIST,AMP_THRESHOLD);
     
            if NumPeaks >1
                sortedLocs = sortedLocs +MINPEAKDIST;
                locs = sortedLocs(1); offset = +4;
                rectPulseCount = 0;
                len = length(sortedLocs);
                for i = 2 : len
                    if (sortedLocs(i) >= (i*locs - offset)) && ...
                        (sortedLocs(i) <= (i*locs + offset))
                       rectPulseCount = rectPulseCount + 1;                    
                    else
                    end
                end

                if rectPulseCount >= floor(len/2)
                   isRectPulse = 1;
                   symratePeakIdx = sortedLocs(1) - MINPEAKDIST;
                end
            else
%                 isRectPulse = 0;
            end                
%             figure; hist(detrendedPartialCM21_dB,100)
        else
        end
        
%         isRectPulse
         %% Is local maximum?
%         [val, idx] = max(detrendedPartialCM21_dB(symratePeakIdx:symratePeakIdx+5));
%         symratePeakIdx = symratePeakIdx + idx - 1;
                
        % edited at 150318
        coarseSymbolRate = partialFreqAxis(symratePeakIdx)

        coarseSamplesPerSymbol  = sam_freq / coarseSymbolRate;
        coarseNumberOfsymbols   = length(rx_sig) / coarseSamplesPerSymbol;   
        
%         serachRange = leftSymRatePeakIdx-2:leftSymRatePeakIdx+2;
%         minVal = min(partialCM21_dB(serachRange));
%         
%         height = maxVal - minVal;
%         if(height > minPeakthreshold) 
%             symratePeak = maxVal;
%             symratePeakIdx =leftSymRatePeakIdx;
%         end
        
%         if strcmp(linearFlag,'nonLinear');   %refer to FMVsDigital.m         
% 
% %             srPower = sum(partialCM21_dB(symratePeakIdx-1 : symratePeakIdx+1));
% %             serachRange = [symratePeakIdx-7 : symratePeakIdx-3 , symratePeakIdx+3: symratePeakIdx+50];
% %             serachRange(serachRange<1) = [];
% %             serachRange(serachRange>length(partialCM21_dB)) =[];
% %             % 'MINPEAKHEIGHT'
% % 
% %             m = mean(partialCM21_dB(serachRange));
% %             s = std(partialCM21_dB(serachRange));
% %             
% %             powerRatio = srPower / m * 100;
% %             threshold =  m + 8*s;
% % 
% %             
% %             anotherPeak = sum( (partialCM21_dB(serachRange)) > minPeakthreshold ); 
% % 
% %             muMinusHill = symratePeak - threshold;
%             if (anotherPeak == 0)
%                 coarseModType = 'FSK';
%             else
%                 coarseModType = 'FM';                     
%                 coarseSamplesPerSymbol  = 0;
%                 coarseNumberOfsymbols   = 0;
%                 coarseSymbolRate        = 0;  
%                 adFlag = 'Analog';
%                 updateLoggingMsg('Modulation type is FM',handles);     
%             end                
%         end
        
        if ~handles.opMode        
            %% update GUI states
            updateLoggingMsg(['# of symbols : ',num2str(coarseNumberOfsymbols)],handles);  
            updateLoggingMsg(['Samples per symbol : ',num2str(coarseSamplesPerSymbol)],handles);
            updateLoggingMsg(['Coarse symbol rate : ',num2str(coarseSymbolRate),'   Mag. of spectrum : ',num2str(symratePeak)],handles);  

            %% plot graph
            cla(handles.axes_monitoring,'reset');
            legend(handles.axes_monitoring,'off');

            plot(partialFreqAxis, partialCM21_dB,'parent',handles.axes_monitoring);
            hold(handles.axes_monitoring,'on');
            plot(partialFreqAxis,spectralTrend,'r','parent',handles.axes_monitoring);
            xlabel('Frequency (Hz)','parent',handles.axes_monitoring);
            ylabel('Magnitude (dB)','parent',handles.axes_monitoring);
            xlim(handles.axes_monitoring,[partialFreqAxis(1), partialFreqAxis(end)]);
            grid(handles.axes_monitoring,'on'); hold(handles.axes_monitoring,'on');        
            plot(partialFreqAxis(symratePeakIdx), ...
            partialCM21_dB(symratePeakIdx),'ro','parent',handles.axes_monitoring);
            legend(handles.axes_monitoring,{'Spectrum','Trend','Symbol rate bin'});
        end
        % legend(handles.axes_monitoring, {'(2,1) cyclic moment','Symbol rate peak'});

    end
    if ~handles.opMode   
        %% update GUI states
        set(handles.pushbutton_coarseToneSpacing,'Enable','on');
        
        set(handles.pushbutton_coarseSymRateEstimation,'BackgroundColor',[222 235 250]/255);
%         set(handles.pushbutton_reSampling,'Enable','on');
        
    end    
    
    %% save data
    handles.sigma_a                = sigma_a;
    handles.coarseSamplesPerSymbol = coarseSamplesPerSymbol;
    handles.coarseNumberOfsymbols  = coarseNumberOfsymbols;
    handles.coarseSymbolRate       = coarseSymbolRate;
    handles.linearFlag             = linearFlag;    % linear digital signal classification bits    
    handles.coarseModType          = coarseModType;
    handles.adFlag                 = adFlag;
%     handles.isFSK                  = isFSK;
%     handles.symStartPoint          = symStartPoint;   % Added at 15.03.05 PM 3:47
    guidata(hObject,handles);
    
% [1].	E. E. Azzouz and A. K. Nandi,  Automatic Modulation Recognition of 
%       Communication Signals. Boston, MA: Kluwer, 1996. 
catch err
    updateLoggingMsg('Error in coarse symbol rate estimation section ',handles);
    rethrow(err);    
end

function handles=pushbutton_reSampling_Callback(hObject, ~, handles)
 updateLoggingMsg('While processing  resampling, please wait ...',handles);
try
    
    adFlag                    = handles.adFlag;
    sam_freq                  = handles.sam_freq;
    rx_sig                    = handles.rx_sig;
    FFTlength                 = handles.FFTlength;
    freqResolution            = handles.freqResolution;
    coarseSamplesPerSymbol    = handles.coarseSamplesPerSymbol;
    coarseNumberOfsymbols     = handles.coarseNumberOfsymbols;
    coarseSymbolRate          = handles.coarseSymbolRate;
    spsLowerBound             = handles.spsLowerBound;
    spsUpperBound             = handles.spsUpperBound;
    previousDecimationFactor  = handles.previousDecimationFactor;
    
    switch adFlag
        case 'NoSignal'
        
        case 'Analog' 
            updateLoggingMsg('While processing decimation, please wait',handles);
            % Decrease sampling rate by integer factor
            p=1;
            q = 512*8 / previousDecimationFactor; % heuristic value                        
            rx_sig = decimate(rx_sig,q);
                        
            sam_freq                  = sam_freq * p / q;
            freqResolution            = freqResolution * p / q;
            currentDecimationFactor   = previousDecimationFactor * q / p;
            
    
        case 'Digital'
             if coarseSamplesPerSymbol < spsLowerBound
                constant = ceil( spsLowerBound / coarseSamplesPerSymbol );            
                upsamplingRatio = 2^nextpow2( constant );   
                p = upsamplingRatio; q=1;
                % Request FPGA to change decimation rate
                % Increase sampling rate by integer factor           
                rx_sig = interp(rx_sig, upsamplingRatio);

            elseif coarseSamplesPerSymbol > spsUpperBound  % what is the proper decimation rate                
                downsamplingRatio = 1;
                while (coarseSamplesPerSymbol / downsamplingRatio) > spsUpperBound  
                    downsamplingRatio = downsamplingRatio + 1;                    
                end
                %constant = floor( coarseSamplesPerSymbol / spsUpperBound	 );
                % 현재의 데시메이션율을 constant보다 큰 2의 배수만큼 늘려야 한다.
                %downsamplingRatio = 2^power; % 2^nextpow2(constant);       
                p = 1;  q=downsamplingRatio;       
                % Decrease sampling rate by integer factor
                %rx_sig = decimate(rx_sig,downsamplingRatio);
                rx_sig = downsample(rx_sig,downsamplingRatio);
             else
                p=1; q=1;
             end
     
            sam_freq                  = sam_freq * p / q;
            freqResolution            = freqResolution * p / q;
            currentDecimationFactor   = previousDecimationFactor * q / p;
            coarseNumberOfsymbols     = length(rx_sig) / coarseSamplesPerSymbol;
            coarseSamplesPerSymbol    = sam_freq / coarseSymbolRate;
           
%             freqAxis = ( -FFTlength/2 : (FFTlength -1)/2 ) / FFTlength * sam_freq;
    
           %% In the multi-rate processing, 
           % remind that the data rate of the signal is constant (unchanged)
            zeroIdx = FFTlength/2+1;
            leftBWIdx = zeroIdx - handles.BWSampleIdx(1);
            rightBWIdx = handles.BWSampleIdx(end)-zeroIdx;            
  
            leftBWIdx = zeroIdx - floor(leftBWIdx * q / p);
            rightBWIdx = floor(rightBWIdx  * q / p) + zeroIdx;
            
            handles.BWSampleIdx = leftBWIdx : rightBWIdx;                          
                       
            updateLoggingMsg([' # of symbols : ',num2str(coarseNumberOfsymbols)],handles);
            updateLoggingMsg(['Samples per symbol : ',num2str(coarseSamplesPerSymbol)],handles);           

        otherwise
          
    end

    viewStartIdx = 1;
    viewEndIdx = length(rx_sig);
    power = nextpow2(currentDecimationFactor);
    
    %% update handles structure 
    handles.rx_sig                          = rx_sig;
    handles.sam_freq                        = sam_freq;
    handles.sam_period                      = 1 / handles.sam_freq;
    handles.freqResolution                  = freqResolution;
    handles.coarseSamplesPerSymbol          = coarseSamplesPerSymbol;
    handles.coarseNumberOfsymbols           = coarseNumberOfsymbols;
    handles.viewEndIdx                      = viewEndIdx;
    handles.viewStartIdx                    = viewStartIdx;    
    handles.previousDecimationFactor        = currentDecimationFactor;
%     handles.freqAxis                        = freqAxis;
    
    guidata(hObject,handles);
    if ~handles.opMode
        %% update GUI states
        set(handles.edit_samFreq,'String',num2str( sam_freq));
        set(handles.edit_freqResolution,'String',num2str(freqResolution));
        set(handles.popupmenu_decimationFactor,'Value',power-4);
        set(handles.edit_inputNumSamples,'String',num2str(length(rx_sig)));
        set(handles.edit_viewStartIdx,'String',num2str(viewStartIdx));
        set(handles.edit_viewEndIdx,'String',num2str(viewEndIdx));

        set(handles.pushbutton_featureExtraction,'Enable','on');

        set(handles.pushbutton_reSampling,'BackgroundColor',[222 235 250]/255);

        updateLoggingMsg('Complete resampling.',handles);    
    end
catch err
    updateLoggingMsg('Error in resampling.',handles);
    rethrow(err);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%   > Classification <   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles=pushbutton_coarseToneSpacing_Callback(hObject, ~, handles)
updateLoggingMsg('While processing coarse tone spacing estimation, please wait ...',handles);
try
    handles=pushbutton_viewSpectrum_Callback(hObject, 0, handles);
    %% Parse the data
    rx_sig            = handles.rx_sig;
    
    coarseBandWidth   = handles.coarseBandWidth;
    coarseSymbolRate  = handles.coarseSymbolRate;
    Pxx               = handles.Pxx;
    
%     coarseModType     = handles.coarseModType; 
    BWSampleIdx       = handles.BWSampleIdx;
    threshold_dB      = handles.BWthreshold;
    
    freqResolution    = handles.freqResolution;
    FFTlength         = handles.FFTlength;

    sam_freq          = handles.sam_freq;
    traceType         = handles.traceType;

    viewStartIdx      = handles.viewStartIdx;
    viewEndIdx        = handles.viewEndIdx;
    linearFlag        = handles.linearFlag;
%     isFSK             = handles.isFSK;
    
    rx_sig = rx_sig(viewStartIdx : viewEndIdx );
    rx_sig = rx_sig(:);

    threshold_dB = threshold_dB-6;
    threshold_lin = db2pow(threshold_dB);
    
    if strcmp(linearFlag,'nonLinear')   
      
        % 옵셋 추정하기 전에 Pxx_dB라서 제대로 안된다.
        % 새로운 스펙트럼 계산이 필요       
        power = 1;
        numTones = 0;
%         figure; plot(real(rx_sig(1:1024))); 
        while numTones <= 1
            if power > 4, break; end
            if power > 1     
                if power == 2 
                    FFTlength = FFTlength / 2;
                end
               %% Calculated power spectral density
                % [Pxx,~] = pwelch(rx_sig,window,noverlap,FFTlength, ...
                %     sam_freq,'twosided');
                window           = hanning( FFTlength );              % Time domain data tapering window
                noverlap         = FFTlength/2;                     % The number of samples for overlap
                freqResolution   = sam_freq / FFTlength;

                numberOfSegments = (length(rx_sig)-noverlap)./(FFTlength-noverlap);
                numberOfSegments = fix(numberOfSegments);

                LminusOverlap = FFTlength-noverlap;
                xStart = 1:LminusOverlap:numberOfSegments*LminusOverlap;
                xEnd   = xStart+FFTlength-1;    
                %accumulatedSxx = zeros(numberOfSegments, FFTlength);
                Sxx = 0; Sxxk=0;
                accumulatedSxx=  zeros(numberOfSegments,FFTlength);    
                switch traceType
                    case 'MaxHold'
                        for i = 1:numberOfSegments     
%                             rx_sig = smooth(rx_sig,FFTlength,'moving');
                            xp = rx_sig(xStart(i):xEnd(i)).^power;                            
                            xw = xp.*window;            % Window the data                                         
                            U = window'*window;         % compensates for the power of the window.
                            Xx = fft(xw,FFTlength);
                            Sxxk = Xx.*conj(Xx)/U;      % Auto spectrum.
                             %Sxx  = Sxx + Sxxk;
                            accumulatedSxx(i,:)  = Sxxk;  % MaxHold Trace in ITU-R SM.443-4
                        end     
                        if size(accumulatedSxx,1) > 1
                            maxHoldSxx = max(accumulatedSxx);        
                        else
                            maxHoldSxx = (accumulatedSxx);
                        end        
                        Pxx = maxHoldSxx;

                    case 'ClearWrite'
                        for i = 1:numberOfSegments            
%                             rx_sig(xStart(i):xEnd(i)) = smooth(rx_sig(xStart(i):xEnd(i)),3,'moving');
                            xp = rx_sig(xStart(i):xEnd(i)).^power;
%                             if i == 1
%                                 figure; plot(real(xp)); 
%                             end
                            xw = xp.*window;                % Window the data
                            U = window'*window;             % compensates for the power of the window.
                            Xx = fft(xw,FFTlength);
                            Sxxk = Xx.*conj(Xx)/U;          % Auto spectrum.
                            Sxx  = Sxx + Sxxk;
                        end     
                        Pxx = Sxx ./ numberOfSegments; % Average the sum of the periodograms
                    otherwise
                end

                Pxx = fftshift(Pxx);          %Pxx_dB =( 10*log10(Pxx) );
                
               %% normalization
               % yMax = max(Pxx_dB); yMin = min(Pxx_dB);
               % Pxx_dB = ((Pxx_dB - yMin) ./ (yMax -yMin)) * 100;
%                 %% calculation BW
%                 num = 1: FFTlength;
% %                 tmp = num( Pxx_dB >= threshold_dB );
%                 tmp = num(Pxx>threshold);
%                 axisLowerIdx = tmp(1);
%                 axisUpperIdx = tmp(end);
% 
%                 BWSampleIdx = axisLowerIdx : axisUpperIdx;
            else
                
            end
            %% calculation BW

            %tmp = num(Pxx>threshold);
            %axisLowerIdx = tmp(1);
            %axisUpperIdx = tmp(end);
            %BWSampleIdx = axisLowerIdx : axisUpperIdx;
            
            %% dB scale
%             Pxx_dB = pow2db(Pxx);
%             logicIdx = Pxx_dB>threshold_dB;
%             partialPxx = Pxx_dB(Pxx_dB>threshold_dB);
%             partialPxx = partialPxx(:); 
%             figure; plot(partialPxx);
                  
           %% linear scale
            logicIdx = Pxx>threshold_lin;
            partialPxx = Pxx(logicIdx);
            partialPxx = partialPxx(:); 

            % for plot
            num = 1: FFTlength;  
            BWSampleIdx = num(logicIdx);
%              %% frequency domain smoothing
%             smoothPartialPxx = smooth(partialPxx, 3,'moving');%length(partialPxx)/8,'moving');
     
            continuousComponents = smooth(partialPxx, floor(length(BWSampleIdx)/4),'moving');%length(partialPxx)/8,'moving');
%             figure; plot(partialPxx); hold on; plot(continuousComponents,'r');
            smoothPartialPxx = partialPxx - continuousComponents;

            %% normalization
%             smoothPartialPxx = partialPxx; %(partialPxx- min(partialPxx)) ./ ...
%                 (max(partialPxx) - min(partialPxx));


            %% peak detection algorithm
            ratioBW2SR = coarseBandWidth / coarseSymbolRate;
            if ratioBW2SR >= 1.8, 
                % Wide tone spacing     
                THRESHOLD = mean(smoothPartialPxx)+2.3*std(smoothPartialPxx); 
            else

                THRESHOLD =  1; %std(smoothPartialPxx); 
            end
            
             %% Remove DC component
%             offset = 7;
            offset = floor(200 / handles.freqResolution  * (power/2)); 
%             offset = floor(  offset/2 );
            len = floor( length(partialPxx) / 2 ); range = len-offset : len+offset;
            smoothPartialPxx(range) = 0; % approximately -330Hz ~ + 330Hz
            
            MINPEAKDIST   = floor(0.8*coarseSymbolRate / freqResolution) ;
            [numTones, sortedToneIdx, threshold]= peakDetection(smoothPartialPxx,MINPEAKDIST,THRESHOLD);

            if numTones == 2
                for i =1 : 2
                     toneMag = smoothPartialPxx(sortedToneIdx(i));

%                     leftRange = sortedToneIdx(i)-floor(MINPEAKDIST/4);
%                     rightRange = sortedToneIdx(i)+floor(MINPEAKDIST/4);
%                     range = leftRange: rightRange;
%                     minVal = min(smoothPartialPxx(range));
%                     toneMag = toneMag - minVal;                
                    if toneMag < 20  %|| localMaxCond > 0
                        numTones = numTones-1;
                    end 
                end
            end

            %% increase power
            power = power*2;
        end

        peakSpacingIdx = diff(sortedToneIdx);
        % mode function returns most frequently occurring value in PeakInterval 
        if isempty(peakSpacingIdx)
            coarseToneSpacing = 0;
        else
            coarseToneSpacing = freqResolution * mean(peakSpacingIdx) / (power/2); 
        end
    	
    
%     		// General information
% 			// 1) 선형변조 신호의 대역폭은 심볼율의 2.5배를 넘지 않는다.
% 			// 2) 비선형변조 신호인 FSK의 대역폭은 심볼율의 8.5배를 넘지 않는다.
% 			//    ex) 심볼율과 대역폭의 비가 가장 큰 프로토콜은 POCSAG으로 6.8배 이다.
% 			// 따라서 심볼율과 대역폭의 비를 계산하면 FSK와 FM을 구분할 수 있다.     
                
        isFM = 0;
        % Check ratio bandwidth to symbol rate
        %if ratioBW2SR > 8,isFM = 1; end

		% Check ratio tone spacing to symbol rate
        ratioTS2SR = coarseToneSpacing / coarseSymbolRate;
        switch(power/2)
            case 1
                if ratioTS2SR < 0.8, isFM = 1; end  
            case 2
                if ratioTS2SR > 0.8, isFM = 1; end  
            case 4
                if ratioTS2SR > 0.4, isFM = 1; end  
            otherwise
        end
        if (coarseToneSpacing <550), isFM=1; end

        % Check the number of Peaks
        if  numTones == 2      %2FSK 
            modulationOrder           = 2;
        elseif 3 <= numTones && numTones <=4  %4FSK
            modulationOrder           = 4;
        elseif 5 <= numTones                 %8FSK  
            modulationOrder           = 8;
        else % if peak spacing is too narrow
            isFM = 1;
        end
      
%         if isFSK>0 % 심볼율 선 스펙트럼의 크기에 따라 FSK 분류 瑩嗤
                 % 톤 스펙트럼 검출에 실패한 경우
                 % from symbol rate est.
%         else            
            if isFM    
                adFlag = 'Analog';
                coarseModType = 'FM';                     
                handles.coarseSamplesPerSymbol  = 0;
                handles.coarseNumberOfsymbols   = 0;
                handles.coarseSymbolRate        = 0;                  
                coarseToneSpacing = 0;
                modulationOrder= 0;                
            else
                adFlag = 'Digital';
                coarseModType = 'FSK';  
            end
%         end
        if ~handles.opMode
            %% plot 
            cla(handles.axes_monitoring,'reset');
            legend(handles.axes_monitoring,'off');
    
            plot(BWSampleIdx,smoothPartialPxx,'parent',handles.axes_monitoring);
            xlim(handles.axes_monitoring,[min(BWSampleIdx), max(BWSampleIdx)]);
            hold(handles.axes_monitoring,'on'); grid(handles.axes_monitoring,'on');
            plot(BWSampleIdx(sortedToneIdx),smoothPartialPxx(sortedToneIdx) ...
                ,'go','parent',handles.axes_monitoring);
            plot(BWSampleIdx,ones(1,length(BWSampleIdx))*threshold ...
                ,'r--','parent',handles.axes_monitoring);  
            xlabel(handles.axes_monitoring,'samples');
            ylabel(handles.axes_monitoring,'magnitude');
            
            updateLoggingMsg('Complete processing ',handles);    
            %% update GUI states
            updateLoggingMsg(['Modulation type is ',coarseModType],handles);
            myMsg = ['Estimated tone spacing : ',num2str(coarseToneSpacing)];
            updateLoggingMsg(myMsg,handles);   
            

        end
    else
        coarseToneSpacing = 0;
        modulationOrder = 0;
        coarseModType = 0;
        adFlag = 'Digital';
        if ~handles.opMode
           %% update GUI states
            updateLoggingMsg('Complete processing ',handles);        
        end
    end
    
    %% update GUI states
    set(handles.pushbutton_featureExtraction,'Enable','on');
    set(handles. pushbutton_coarseToneSpacing,'BackgroundColor',[222 235 250]/255);
    
    %% save handle structure
    handles.coarseToneSpacing       = coarseToneSpacing;
    handles.modulationOrder         = modulationOrder;
    handles.coarseModType           = coarseModType;
    handles.adFlag                  = adFlag;
    guidata(hObject,handles);
catch err
    updateLoggingMsg('Error in coarse tone spacing estimation section ',handles);
    rethrow(err);    
end

function handles=pushbutton_featureExtraction_Callback(hObject, ~, handles)

updateLoggingMsg('While processing feature extraction, please wait ...', handles);
try

    viewStartIdx    = handles.viewStartIdx;
    viewEndIdx      = handles.viewEndIdx;
    limitation      = handles.limitation;
    rx_sig          = handles.rx_sig;
    sigma_a         = handles.sigma_a;
    rx_sig = rx_sig(viewStartIdx:viewEndIdx);

    %% Normalization
    rx_sig = sqrt(length(rx_sig)) * rx_sig / norm(rx_sig,2);  
    %% Remove mean
    rx_sig = rx_sig - mean(rx_sig);

    updateLoggingMsg('Loading reference data',handles);

    samplesPerFrame = 32768;  %samplesPerFrame=16384;
    %% Calculate the number of segments
    if length(rx_sig) < samplesPerFrame
        samplesPerFrame = length(rx_sig); numSegment=1;
    else
        numSegment = floor( length(rx_sig) / samplesPerFrame );             
    end

    rmsDifference         = zeros(numSegment,1);   % For classification between ASK and PSK,QAM
    cc20             = zeros(numSegment,1);
    cc40             = zeros(numSegment,1);
    updateLoggingMsg('Extracting features from the input IQ data',handles);
    segIndex = 1: samplesPerFrame;

    for frameCount=1:numSegment

        % divide signal with size of samplesPerFrame
        partialSig = rx_sig(segIndex); 
        segIndex = segIndex+samplesPerFrame;        

        %% Instantaneous information
        % instantaneous amplitude
        a_amplitude = abs(partialSig);            
        % instantaneous phase + 2*pi*fc*t
        phi_phase = angle(partialSig);  % between +- pi

        c_k(1)=0; 
        % make unwrapped angle
        for i=2:samplesPerFrame
            if phi_phase(i)-phi_phase(i-1)>pi
                c_k(i)=c_k(i-1)-2*pi;
            elseif phi_phase(i)-phi_phase(i-1)<-pi
                c_k(i)=c_k(i-1)+2*pi;
            else c_k(i)=c_k(i-1);
            end
        end

        phi_phase = phi_phase(:);
        c_k = c_k(:);
        phi_uw_phase=phi_phase+c_k;
        % subtract the linear phase component from unwrapped angle
        % instantaneous phase, remove the carrier component from analytic signal.
        %  phi_NL_phase=phi_uw_phase-2*pi*Param.carrier_freq*t_vec;
        phi_NL_phase=phi_uw_phase;        
        % equal code        
        % a_phase =angle( analyticSignal.*exp(-1i*2*pi*Param.carrier_freq*t_vec));
        y = phi_NL_phase ;
        x = 1 : samplesPerFrame;
        v =[x.' , ones(samplesPerFrame,1)];
        % increasing order
        p = v\y;    % refer to iLinearFit.m 
%             figure; plot(y(1:3000));
        % Same as 
        % [Q,R] = qr(v,0);
        % p2 = full(R\(Q'*y));    % Same as p = D*A\(D*y);
        haty = v*p;
        rmsDifference(frameCount) =  sqrt(mean((y - haty).^2) );  

       %% high order statistics
        cc20(frameCount) = cm(partialSig,2,0,handles);
        cc40(frameCount) = cm(partialSig,4,0,handles);
        cc40(frameCount) = cc40(frameCount) - 3* cc20(frameCount)^2;
    end 
    updateLoggingMsg('Outliers',handles);

    %% except outlier for inputData  added 141218        
    rmsDifference = dataConditioning(rmsDifference);
    cc20 = dataConditioning(cc20); % for BPSK vs other linear modulation
    cc40 = dataConditioning(cc40);
    [~,sigma_a_middle] =  dataConditioning(sigma_a);
    %% Save data
    handles.limitation = limitation;
    % handles.featureLabel = featureLabel;

    handles.cm40                  = cc40;
    handles.cc20                  = cc20;         % for BPSK vs other linear modulation
    handles.rmsDifference         = rmsDifference;
    handles.numSegment            = numSegment;
    updateLoggingMsg(['1.|CC20| : ',num2str(median(cc20)), ...
        '   2.|CC40| : ',num2str(median(cc40)), ...
        '   3.std. of inst. amp. : ',num2str((sigma_a_middle)) ],handles);        
    guidata(hObject,handles);
    if ~handles.opMode
        %% plot
%             plotRefDataNfeatureSet(handles);
        cla(handles.axes_monitoring,'reset');
        legend(handles.axes_monitoring,'off');
        
        fSet=[cc20, cc40, rmsDifference/1000];
        bar(handles.axes_monitoring,median(fSet));
%         set(handles.axes_monitoring,'Xtic',['cc20','cc40','rms'])
        %% update GUI states
        updateLoggingMsg('Complete feature extraction',handles);
        set(handles.pushbutton_fusionCenter,'Enable','on');

        set(handles.pushbutton_featureExtraction,'BackgroundColor',[222 235 250]/255);
        set(handles.uipushtool_viewForwardFig,'Enable','on');
        set(handles.uipushtool_viewBackFig,'Enable','on');

        % handles.featureWeights = eval(strcat('[',get(hObject,'String'),']'));
        % set(handles.edit_bandwidthConstant,'String',num2str(featureWeights.'));
        guidata(hObject,handles);
    end

catch err
    updateLoggingMsg('Error in feature extraction',handles);
    rethrow(err);    
end

function handles = pushbutton_fusionCenter_Callback(hObject, ~, handles)

updateLoggingMsg('While processing recognition, please wait ...',handles);

Lmod            = handles.Lmod;
adFlag          = handles.adFlag;
linearFlag      = handles.linearFlag;
% coarseModType   = handles.coarseModType;
% modulationOrder = handles.modulationOrder;
% rmsDifference        = handles.rmsDifference;
% cc20            = handles.cc20;
cm40            = handles.cm40;
sigma_a         = handles.sigma_a;
isDehoppingSig  = handles.isDehoppingSig;

if isDehoppingSig == 1
    % We assume modulation types of FH sig always is FM 
    selectedLmod = Lmod(1:2); 
    modulationOrder = 0;
    totalPercents = zeros(1,2);
    totalPercents( 2 ) = 100;    
else
    switch adFlag
        case 'NoSignal'

        case 'Analog' 
            % class label        
            selectedLmod = Lmod(1:2); 
            coarseModType = handles.coarseModType;
            switch coarseModType
                case 'AM'  %AM                            
                    modulationOrder = 0;
                    totalPercents = zeros(1,length(selectedLmod));
                    totalPercents( 1 ) = 100;
                case 'FM' %FM                        
                    modulationOrder = 0;
                    totalPercents = zeros(1,length(selectedLmod));
                    totalPercents( 2 ) = 100;
            end    
        case {'Digital'}
            switch linearFlag
                case 'Linear'

                    selectedLmod = Lmod([3 7:length(Lmod) ]);      
                    modulationOrderTable = [2,2,4,8,32];
                    if median(handles.rmsDifference) < 1.5 % refer to testLMScurve.m
                        coarseModType = selectedLmod{1};     % ASK
                        modulationOrder = 2;
                        totalPercents = zeros(1,length(selectedLmod));
                        totalPercents( 1 ) = 100;
                    else
                        if median(handles.cc20) > 0.5                        
                            coarseModType   = selectedLmod{2};   % BPSK
                            modulationOrder = 2;                 
                            totalPercents = zeros(1,length(selectedLmod));
                            totalPercents( 2 ) = 100;
                        else
                            if sigma_a > 0.37 % refer to SigmaForLinearMod.fig at 150420
                                % QAM 
                                modCount =[0 0 0 0 1];
                            else
                                if cm40 > 0.14 % refer to CC40 forMPSK.fig at 150420
                                   % QPSK 
                                   modCount =[0 0 1 0 0];
                                else
                                   % 8PSK
                                   modCount =[0 0 0 1 0];
                                end                            
                            end
                            % counting
    %                         x = modNumber(oneColumnFeatureIdx) ;
                            % index : 1 : 2ASK, 2:BPSK, 3:QPSK, 4:8PSK, 5:QAM
    %                         [modCount, ~] = hist(x,1:length(selectedLmod));                   

                            % Get percentage 
                            total = sum(modCount);
                            totalPercents = 100*modCount./total; 

                            % Get maximum similarity
                            [~, maxIdx]     = max(totalPercents);
                            coarseModType   = selectedLmod{maxIdx};     
                            modulationOrder = modulationOrderTable(maxIdx);
                        end
                    end       
                case 'nonLinear'
                    coarseModType   = handles.coarseModType;
                    modulationOrder = handles.modulationOrder;
                    selectedLmod    = Lmod([4,5,6]);                                                   
                    totalPercents   = zeros(1,length(selectedLmod));
                    totalPercents( log2(modulationOrder) ) = 100;
                otherwise
            end                
       otherwise
    end  % for adFlag      
end % for isDehoppingSig?
%% plot grpah
cla(handles.axes_monitoring,'reset');
legend(handles.axes_monitoring,'off');

bar(handles.axes_monitoring, 1:length(selectedLmod), totalPercents);
xlim(handles.axes_monitoring,[ 1-0.5 length(selectedLmod)+0.5]);
ylabel(handles.axes_monitoring,'Similarity (%)');
set(handles.axes_monitoring,'XTickLabel',selectedLmod);
%%%%%%%%%%%%%%%%%%%%%%%%%%% > result of PCA < %%%%%%%%%%%%%%%%%%%%%%%%%% 
% trnData     = refDataPrincipleComponent(modSigRange,:,estimatedSNRIdx);
% [featureIdx, distanceFromFeature] = knnsearch(trnData(:,1:3),newPointsPCA,'k',kNeighbors); 
% oneColumnFeatureIdx = featureIdx(:); 

 %% update GUI states  
set(handles.edit_viewBW,'String',num2str(floor(handles.coarseBandWidth)));
set(handles.edit_viewSNR,'String',num2str(floor(handles.estimatedSNR)));
set(handles.edit_viewCFO,'String',num2str(floor(handles.coarseFreqOffset)));
if strcmp(adFlag,'Digital')
    set(handles.edit_viewSymRate,'String',num2str(floor(handles.coarseSymbolRate)));
     set(handles.edit_viewNumSym,'String',num2str(floor(handles.coarseNumberOfsymbols)));
    set(handles.edit_viewSPS,'String',num2str(floor(handles.coarseSamplesPerSymbol)));
   % if strcmp(linearFlag,'nonLinear')
      set(handles.edit_viewToneSpacing,'String',num2str(floor(handles.coarseToneSpacing)));
  % end
end
updateLoggingMsg('Complete recognition',handles);
set(handles.pushbutton_fusionCenter,'BackgroundColor',[222 235 250]/255);

%% save data
handles.coarseModType   = coarseModType;
handles.modulationOrder = modulationOrder;
guidata(hObject,handles);

%% Enable demodulation process
set(handles.pushbutton_demod,'Enable','on');
%%%%%%%%%%%%%%%%%%%%%%%%%%   <Classification>   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Spectral-line methoud based on fourth-power nonlinearity 
    % NthOrder = 8;
    % FFTlength = 2^(nextpow2(length(rx_sig)));
    % NthPoweredSig = rx_sig.^NthOrder;
    % spectrum =  (fft(NthPoweredSig,FFTlength)) ;
    % abs_spectrum = abs(spectrum);
    % angle_spectrum = angle(spectrum);
    % abs_spectrum_dB = 10*log10(abs_spectrum);
    % 
    % [val, indexFreqOffset]=max(abs_spectrum_dB);
    % 
    % % the N-th powered spectrum in the left side from center moves to center
    % % Also, the spectrum in the right side from center moves to center
    % % Therefore, we need below condition.  date : 141105
    % 
    % ttIdx = 0 : length(rx_sig)-1;
    % rx_sig = rx_sig(:);
    % ttIdx = ttIdx(:);
    % freqAxis = (0 : FFTlength-1) / FFTlength;
    % 
    % if indexFreqOffset > FFTlength /2
    %     indexFreqOffset = FFTlength - indexFreqOffset;
    %     carrierFreqOffset = freqAxis(indexFreqOffset) / NthOrder;
    %     phaseOffset = angle_spectrum(indexFreqOffset) / NthOrder;
    %     rx_sig = rx_sig.*exp(1i*(2*pi*carrierFreqOffset*ttIdx + phaseOffset) );
    % else
    %     carrierFreqOffset = freqAxis(indexFreqOffset) / NthOrder;
    %     phaseOffset = angle_spectrum(indexFreqOffset) / NthOrder;
    %     rx_sig = rx_sig.*exp(-1i*(2*pi*carrierFreqOffset*ttIdx + phaseOffset) );
    % end

    %% Assume that sampling frequency is 1
    % freqOffsetRemoved_abs_spectrumd_dB = ...
    %     10*log10( abs(fft(rx_sig.^NthOrder,FFTlength)) );
    
% function pushbutton_fineCarrierOffsetEst_Callback(hObject, eventdata, handles)
% 
% updateLoggingMsg('While processing fine carrier offset estimation, please wait ...',handles);
% 
% coarseModType =handles.coarseModType;
% rx_sig = handles.rx_sig ;
% % estBW = handles.coarseBandWidth;
% sam_freq = handles.sam_freq;
% viewEndIdx   = handles.viewEndIdx;
% viewStartIdx = handles.viewStartIdx;
% 
% rx_sig = rx_sig(viewStartIdx : viewEndIdx);
%         
% switch coarseModType
%     case {'AM','FM'}
%         estimatedFineCarrierOffset = 0;
%         updateLoggingMsg('Skip',handles);
%     case {'2FSK','4FSK','8FSK'}
%         estimatedFineCarrierOffset = 0;
%         updateLoggingMsg('Skip',handles);      
%     otherwise
%       % Step 2 : Estimate (2,1) cyclic moment       
%         [normMag, resolution, NFFT] = calculateCyclicMoment(rx_sig,2,1,handles);   
%         normMag = fftshift(normMag);
% 
%         halfNormMag = normMag(1:end/2);
%         % xcorrNormMag = xcorr(normMag); % too many time consumption
%         freqAxis = (0:NFFT/2-1) /NFFT * sam_freq ;
% %         freqAxis = (-NFFT/2:NFFT/2-1) /NFFT * sam_freq ;
%         % dcFreeHalfNormalizedMag = normalizedMag(dcOffset:end/2);
%         hlafNormMag_dB = 10*log10(halfNormMag);
%         
%         % Step 5 : Search the global maximum peak within range        
%         [offsetPeakVal, offsetPeakIdx] = max(hlafNormMag_dB);
%         estimatedFineCarrierOffset = freqAxis(offsetPeakIdx);
%         
%         % plot graph
%         cla(handles.axes_monitoring,'reset');
%         legend(handles.axes_monitoring,'off');
%         
%         plot(freqAxis, hlafNormMag_dB, 'parent',handles.axes_monitoring);          
%         hold(handles.axes_monitoring,'on');
%         plot(freqAxis(offsetPeakIdx), normMag(offsetPeakIdx),...
%             'ro','parent',handles.axes_monitoring);        
%         legend(handles.axes_monitoring, {'(2,1) cyclic moment','Carrier offset peak'});
%         
% end
% 
% %% save data
% handles.estimatedFineCarrierOffset = estimatedFineCarrierOffset;
% guidata(hObject,handles);
% 
% %% update GUI states
% set(handles.pushbutton_fineCarrierOffsetEst,'BackgroundColor',[222 235 250]/255);
% updateLoggingMsg(['Estimated fine carrier offset is ',num2str(estimatedFineCarrierOffset), ' Hz'],handles);

function [normalizedPartialData,normalizedSigFeature,minMaxData] = normalization(refData, partialSig,num)
%% normzlize the feature set
[dataRow, dataCol, dataZ] = size(refData);

[sigRow, sigCol] = size(partialSig);

% normalization data range between 0 to 1
normalizedPartialData = zeros(dataRow, dataCol, dataZ);
minData = zeros(1, dataZ); maxData = zeros(1, dataZ);
normalizedSigFeature = zeros(sigRow,dataZ);
for zloop = 1 : dataZ   % SNR loop
    minData(zloop) = min(min(refData(:,:,zloop)));    maxData(zloop) = max(max(refData(:,:,zloop)));
    normalizedPartialData(:,:,zloop) = (refData(:,:,zloop) - minData(zloop)) ...
        ./ (maxData(zloop) - minData(zloop));
 
    normalizedSigFeature(:,zloop) =  (partialSig- minData(zloop)) / ...
        (maxData(zloop) - minData(zloop));
end

tmp = [minData;maxData];
minMaxData = tmp(:);


function [output, q2] = dataConditioning(input)
% source from 
% http://alex.bikfalvi.com/research/advanced_matlab_boxplot/

[row,col,zz] = size(input);
output = zeros(row,col,zz);

for zloop = 1 : zz
    for jloop = 1 : col
        oneColumnData = input(:, jloop, zloop);
        % rank the data
        [y, yIdx] = sort(oneColumnData);

        % compute 50th percentile (second quartile)
        q2 = nanmedian(y);

        % compute 25th percentile (first quartile)
        q1 = nanmedian(y(y<=q2));

        % compute 75th percentile (third quartile)
        q3 = nanmedian(y(y>=q2));

        % compute Interquartile Range (IQR)
        IQR = q3-q1;

%         fl = min(y(y>=q1-1.5*IQR)); % Upper Adjacent
%         fu = max(y(y<=q3+1.5*IQR)); % Lower Adjacent

        % 1.5 corresponds to approximately +/-2.7sigma and 99.3 coverage if
        % the data are normally distributed.
%         ol = y(y<q1-1.5*IQR);       % Lower outlier
%         ou = y(y>q3+1.5*IQR);       % Upper outlier

        oneColumnData(oneColumnData<q1-1.5*IQR) = q2;
        oneColumnData(oneColumnData>q3+1.5*IQR) = q2;

        output(:, jloop,zloop) = oneColumnData;
    end 
end

% We cannot use below code in that distribution of data is unknown
%
% get rid of samples outside the range of 95% confidential interval  
% dataL = mean_Feature - 1.96*std_Feature;
% dataR = mean_Feature + 1.96*std_Feature;
% condition = repmat(dataL,limitation, 1)<= partialData ...
%     & partialData <= repmat(dataR,limitation,1);
% partialData(~condition) =0;
% %%%%%         figure; hist(nonzeros(nandi_gamma_max(1:limitation,1,estimatedSNR)),50)
% outData = partialData;

function output = momNcum(rx)
% Description : compute the set of cumulants
%
% edited on Saturday June 2013 by AWH

% remove the mean value from input signal
rx = rx-mean(rx);
m20 = M(rx,2,0);

% signal power
% m21 = M(rx,2,1); 

m40 = M(rx,4,0);
% c40 = m40 - 3*m20.^2;

% c41 = M(x,4,1) - 3*M(x,2,0).*M(x,2,1);

% c42 = M(rx,4,2) - abs( m20 ).^2 - 2*m21.^2;
% c42n = c42./ ( (c21-c21g)^2 );

c60 = M(rx,6,0) - 15*m20.*m40 + 30*m20.^3;

% m41 = M(rx,4,1);
% c61 = M(rx,6,1) - 5*m21.*m40 - 10*m20.*m41 + ...
%       30*m20.^2*m21;    

% c62 = M(x,6,2) - 6*M(x,2,0).*M(x,4,2) - 8*M(x,2,1).*M(x,4,1) - ...
%     M(x,2,2).*M(x,4,0) + 6*M(x,2,0).^3 + ...
%     24*M(x,2,1).^2 * M(x,2,0);

% [1]
% c63 = M(rx,6,3) - 9*m21.*M(rx,4,2) + 12*m21.^3 - ...
%       3*m20.*M(rx,4,3) -3*M(rx,2,2).*m41 + ...
%       18*m20*m21*M(rx,2,2);
% equvalant to 
%  [2] c63 = M(x,6,3) - 9*M(x,2,1).*M(x,4,2) + 12*M(x,2,1).^3 + ...
%     12*abs(M(x,2,0)).^2.*M(x,2,1) 
%
%  [3] cc63 = M(x, 6, 3 ) - 6*M(x, 2, 2).*M(x, 4, 3)...
%     - 9*M(x, 2, 1).*M(x, 4, 2) + 18*M(x, 2, 2).^2.*M(x, 2,1 ) ...
%     + 12*M(x, 2, 1).^3; 

% c80 = M(rx,8,0) - 28*M(rx,6,0).*m20 - 35*m40.^2 + ...
%       420*m40.*m20.^2 -630*m20.^4;
% c80= M(x,8,0) - 35*M(x,4,0).^2 + ...
%    420*M(x,4,0).*M(x,2,0).^2 -630*M(x,2,0).^4;  % [4]

% c81 = M(rx,8,1)-35*m40.*m41-630*m20.^3.*m21 + ... 
%       210*m40.*m20.*m21  + 210*m20.*m41;
% 
% c82 = M(rx,8,2)-15*m40.*M(rx,4,2)-20*m41.^2 + ... %[5]
%     30*m40.*m20.^2 + 60*m40.*m21.^2 + ...
%     240*m41.*m21.*m20 + 90*M(rx,4,2).*m20.^2 - ...
%     90*m20.^4 - 540*m20.^2*m21.^2;
% 
% c83 = M(rx,8,3)-5*m40.*m41-30*M(rx,4,3).*M(rx,4,2) + ...  %[5]
%     90*M(rx,4,3).*m20.^2 + 120*M(rx,4,3).*m21.^2 + ...
%     180*M(rx,4,2).*m21.*m20 + 30*m40.*m20.*m21 - ...
%     270*m20.^3.*m21 - 360*m21.^3*m20;
% 
% c84 = M(rx,8,4)-M(rx,4,4).^2-18*M(rx,4,2).^2-16*M(rx,4,3).^2 -...  %[4]
%     54*M(rx,2,2).^4 - 144*m21.^4-432*M(rx,2,2).^2.*m21.^2 + ...
%     12*M(rx,4,4).*M(rx,2,2).^2 + 96*M(rx,4,3).*m21.*M(rx,2,2) + ...
%     144*M(rx,4,2).*m21.^2 + 72*M(rx,4,2).*M(rx,2,2).^2 + ...
%     96*M(rx,4,3).*M(rx,2,2).*m21; 

% c85 = M(rx,8,5)-5*M(rx,4,4)*M(rx,4,3)-30*M(rx,4,3)*M(rx,4,2) + ...
%     90*M(rx,4,3)*M(rx,2,2).^2 + 120*M(rx,4,3)*m21.^2 + ...
%     180*M(rx,4,2)*m21*M(rx,2,2) + 30*M(rx,4,4)*M(rx,2,2)*m21 - ...
%     270*M(rx,2,2).^3*m21 - 360*m21.^3*M(rx,2,2);
% 
% c86 = M(rx,8,6)-15*M(rx,4,4).*M(rx,4,2)-20*M(rx,4,3).^2 + ... %[4]
%     30*M(rx,4,4).*M(rx,2,2).^2 + 60*M(rx,4,4).*m21.^2 + ...
%     240*M(rx,4,3).*m21.*M(rx,2,2) + 90*M(rx,4,2).*M(rx,2,2).^2 - ...
%     90*M(rx,2,2).^4 - 540*M(rx,2,2).^2.*m21.^2;
% 
% c87 = M(rx,8,7) -35*M(rx,4,4).*M(rx,4,3)-630*M(rx,2,2).^3*m21+...
%     210*M(rx,4,4)*M(rx,2,2)*m21 + 210*M(rx,2,2)*M(rx,4,3);
% 
% c88 = M(rx,8,8) - 35*M(rx,4,4).^2 + ...
%         420*M(rx,4,4).*M(rx,2,2).^2 -630*M(rx,2,2).^4;  %[4]
        
% [1] "Automatic Modulation Classification Using Combination of
%     Genetic Programming and KNN", 2012
% [2] "Automatic Modulation Classification Sixth-order Cumulant Features 
%      as a Solution for Real-world Challenges", 2012
% [3] "Performance Analysis and Optimization of Novel High-Order Statistic 
%      Features in Modulation Classification", 2008
% [4] "Classification of Digital Modulation Types in Multipath
%      Environments", Naval Postgraduate School, 2001
% output = abs( [c40 c42 c61 c63 c80 c84] );
output = abs( c60 );

function m=M(y,p,q)
m=mean(y.^(p-q) .* conj(y).^q );

function res=cm(x,p,q,Param) % Cyclic moment function

antiCycleAliasing(Param, p) % [2]

L=length(x);              % signal length
if L > 65536, L = 65536; end

nonlinearTransform = x.^(p-q) .* conj(x).^q ;

% Xilinx FFT IP core supports the maximum 2^16 point FFT
NFFT = 2^nextpow2(L);        % FFT point 

% if NFFT > cmNFFT, NFFT = cmNFFT; end;

Xk = 1/L * fft( nonlinearTransform, NFFT )  ;  

% halfXk = Xk(1:NFFT/2);
% halfXk = Xk;
% By choosing the CC magnitude, the rotation caused by a fixed phase shift
% and signal delay disappears [3].
MagXk = abs(Xk); %if p == 2 , figure; plot(MagXk);end
[res, maxIdx] = max( MagXk );      % easy way, it should be improved

function fig_AMRtestbed_WindowButtonDownFcn(hObject, eventdata, handles)
axes_area = get(handles.axes_monitoring,'Position');
try
    MouseLoc =  get(handles.fig_AMRtestbed,'currentpoint');
    
    axesLeftSide = axes_area(1) + axes_area(3);
    axesRightSide = axes_area(1);
    axesTop =  axes_area(2) + axes_area(4);
    axesBottom =  axes_area(2);
    
    
    if (MouseLoc(1) < axesLeftSide ) &&...
        (MouseLoc(1) > axesRightSide)  &&...
        (MouseLoc(2) < axesTop) &&...
        (MouseLoc(2) > axesBottom)
    
        axesMouseLoc = get(handles.axes_monitoring,'currentpoint');
        currentXpointOnAxes = axesMouseLoc(1,1);
        currentYpointOnAxes = axesMouseLoc(1,2);
        updateLoggingMsg( ...
           ['X-point:',num2str(currentXpointOnAxes),'    Y-point(xdB point):',num2str(currentYpointOnAxes)],handles);

        handles.currentXpointOnAxes = currentXpointOnAxes;
        handles.currentYpointOnAxes = currentYpointOnAxes;        
        
    else
    end
    guidata(hObject,handles);
catch
    updateLoggingMsg('Error in fig_AMRtestbed_WindowButtonDownFcn',handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot section
function plotGraph(handles)
% clear previous graph    
cla(handles.data.reqAxes);
legend(handles.data.reqAxes,'off');

x = handles.data.x; x=x(:).';
y = handles.data.y; y=y(:).';       
switch handles.units
    case 'dBm' % dBm
        plot(x, 10*log10(y) + 30 ,'parent',handles.data.reqAxes);
        xlabel('Frequency (Hz)','parent',handles.data.reqAxes)
        ylabel('Magnitude (dBm/Hz)','parent',handles.data.reqAxes)
        xlim(handles.data.reqAxes,[x(1), x(end)]);
        grid(handles.data.reqAxes,'on'); 

    case 'dB' % dB
        plot(x, 10*log10(y),'parent',handles.data.reqAxes);
        xlabel('Frequency (Hz)','parent',handles.data.reqAxes);
        ylabel('Magnitude (dB/Hz)','parent',handles.data.reqAxes);
        xlim(handles.data.reqAxes,[x(1), x(end)]);
        grid(handles.data.reqAxes,'on'); 

    case 'Watts' % Watts
        plot(x, y,'parent',handles.data.reqAxes);
        xlabel('Frequency (Hz)','parent',handles.data.reqAxes)
        ylabel('Magnitude (Watts/Hz)','parent',handles.data.reqAxes)
        grid(handles.data.reqAxes,'on'); 
        xlim(handles.data.reqAxes,[x(1), x(end)]);
    otherwise
end

%% update monitoring message box
function updateLoggingMsg(myMsg,handles)
oldmsgs = cellstr(get(handles.text_monitor_msg,'String'));
if length(oldmsgs) >=5 
    set(handles.text_monitor_msg,'String',[...
        oldmsgs{end-3};oldmsgs{end-2};oldmsgs{end-1};oldmsgs{end};{myMsg}] );
else
    set(handles.text_monitor_msg,'String',[oldmsgs;{myMsg}] );
end
pause(0.000001);         

function handles=plotRefDataNfeatureSet(handles)

% currentPlot           = handles.currentPlot;
% axes_monitoring       = handles.axes_monitoring;
% estimatedSNRIdx       = handles.estimatedSNRIdx;
% normalizedFeaturesSet = handles.normalizedFeaturesSet;
% limitation            = handles.limitation;
% modIdx                = handles.refModIdx;
% featureLabel          = handles.Fset;
% Lmod                  = handles.numOfLienarModScheme;
% newPointsPCA          = handles.newPointsPCA;
% featureLabel          = handles.featureLabel;
% numSegment            = handles.numSegment;
% refData = handles.refData;
% refDataPrincipleComponent = handles.refDataPrincipleComponent;
% projectionMat = handles.projectionMat;

% markerSet      = {'o','+','*','.','s','d','^','v','>','p','h','<','o','+'};
% modType = handles.modType;

%% load reference data
% switch refDataType
%     case 'MATLAB'
%         load refData_nandi.mat;      % 20141006
%         load refData_hosC.mat; 
%         load refData_hosCC.mat; 
%         load refData_circ.mat;      % 20141006

%     case 'E4438C'
% %                   load refData_FromGenerator_150131.mat
% %                   load refData_FromGenerator_150203.mat % 2ASK 제거
%           load refData_FromGenerator_150407.mat % cc20 추가, gammax_max, circular statistics 제거
%     otherwise
% end
        
% cla(axes_monitoring,'reset');
% legend(axes_monitoring,'off');
% set(axes_monitoring,'Position',[10.2, 9.0, 102.8, 33.076923 ]);
% hold(axes_monitoring,'on');
% xx = repmat(1:size(normalizedFeaturesSet,2),size(normalizedFeaturesSet,1),1);
% xx = repmat(1:12,size(normalizedFeaturesSet,1),1);
% xx2 = repmat(modIdx,limitation,1);
% 
% featureName ={'nandi_sigma_a','hoc_cc20',...
%     'hoc_cc40'};

%% plot reference data
  
% str=['handles.',featureName{currentPlot},'(1:limitation,modIdx,estimatedSNRIdx)'];
%    boxplot(handles.hoc_cc40(1:limitation,modIdx,estimatedSNRIdx),Lmod,...
%         'labelorientation','horizontal','parent',axes_monitoring);  
%     plot(xx,repmat(normalizedFeaturesSet(:,estimatedSNRIdx,6),1,12), 'g.','parent',axes_monitoring); 
%     ylabel(featureLabel{6});   

% plot(xx2,eval(str),'.','parent',axes_monitoring);    
% plot(modIdx,repmat(normalizedFeaturesSet(:,estimatedSNRIdx,currentPlot),1,length(modIdx)), 'go','parent',axes_monitoring);
% ylabel(featureLabel{currentPlot}); xlabel({'QPSK=1, 8PSK=2, QAM=3'});


% grid(axes_monitoring,'on');
% set position of axes
% set(axes_monitoring,'Position',[10, 13.15,100.4,34.23]);

function uipushtool_viewForwardFig_ClickedCallback(hObject, eventdata, handles)
% currentPlot = handles.currentPlot;
% Fset = handles.Fset;
% if currentPlot == length(Fset)
%     currentPlot = 0;
% else
%     currentPlot = currentPlot +1;
% end
% handles.currentPlot = currentPlot;
% guidata(hObject,handles);

%% plot
plotRefDataNfeatureSet(handles);

function uipushtool_viewBackFig_ClickedCallback(hObject, eventdata, handles)
% currentPlot = handles.currentPlot;
% Fset = handles.Fset;
% if currentPlot == 0
%     currentPlot = length(Fset);
% else
%     currentPlot = currentPlot - 1;
% end
% 
% handles.currentPlot = currentPlot;
% guidata(hObject,handles);

%% plot
plotRefDataNfeatureSet(handles);

function pushbutton_demod_Callback(hObject, eventdata, handles)
AMRtestbedHandles.rx_sig                   = handles.rx_sig;
AMRtestbedHandles.coarseModType            = handles.coarseModType;
AMRtestbedHandles.sam_freq                 = handles.sam_freq;                   % sampling frequency
AMRtestbedHandles.coarseBandWidth          = handles.coarseBandWidth;
AMRtestbedHandles.coarseFreqOffset         = handles.coarseFreqOffset;
AMRtestbedHandles.coarseSamplesPerSymbol   = handles.coarseSamplesPerSymbol;
AMRtestbedHandles.coarseSymbolRate         = handles.coarseSymbolRate;
AMRtestbedHandles.coarseToneSpacing        = handles.coarseToneSpacing;
AMRtestbedHandles.modulationOrder          = handles.modulationOrder;
AMRtestbedHandles.previousDecimationFactor = handles.previousDecimationFactor;
AMRtestbedHandles.coarseModType            = handles.coarseModType;
save('AMRtestbedHandles', 'AMRtestbedHandles');
Demodulation();

% function pushbutton_peakDetection_Callback(hObject, eventdata, handles)
% % hist-based peak detection
% Pxx                 = handles.Pxx;
% NMag                = Pxx./max(Pxx); %normalized magnitude
% Param.MinimumHeight = 1;
% [NumPeaks, estimatedPeakDistance] = peakDetection_140811(NMag.',100,Param);
% 
% if NumPeaks == 1 
%    set(handles.text_2ASK,'FontWeight','bold','ForegroundColor',[1 0 0 ]); 
% else
%    set(handles.text_1stPDTo2ndPD,'FontWeight','bold','ForegroundColor',[1 0 0 ]); 
% end
% 
% set(handles.pushbutton_peakDetection,'BackgroundColor',[222 235 250]/255);
% % save data
% guidata(hObject,handles);

% figure;
% subplot(2,1,1), plot(percentage,sorted_smoothedPxx);
% ylabel('Pxx (dB)'); xlabel('x%'); grid on; hold on;
% subplot(2,1,2),plot(percentage(beginIdx:end),dist,'r',percentage(beginIdx+minIdx),dist(minIdx),'go'); 
% ylabel('Distance from line'); xlabel('Percentage rate of between integrated power and whole power(%)');
% title(['SNR =',num2str(snr_dB(estimatedSNR))]);


% function handles=pushbutton_apen_Callback(hObject, eventdata, handles)
% noisy_rx_sig    = handles.noisy_rx_sig;
% lineSet         = handles.lineSet;
% Nmod            = handles.Nmod;
% samplesPerFrame = handles.samplesPerSegment;
% snr_dB          = handles.featureSNR_dB ;
% numSegment      = floor( length(noisy_rx_sig) / samplesPerFrame );
% fSet            = {'Real ApEn','Imaginary ApEn'};
% samplesPerFrame = handles.samplesPerSegment ; % user optional
% modIdx          = handles.refModIdx;
% limitation      = handles.limitation;
% 
% if length(noisy_rx_sig) < samplesPerFrame
%     samplesPerFrame = length(noisy_rx_sig); numSegment=1;
% else
%     numSegment = floor( length(noisy_rx_sig) / samplesPerFrame ); 
% end
% 
% realApEnData = zeros( numSegment,1);
% imagApEnData = zeros( numSegment,1);
% 
% for frameCount=1:numSegment
%     % divide signal with size of nSegment
%     segIndex = 1+(frameCount-1)*samplesPerFrame : frameCount*samplesPerFrame;  
%     partialData = noisy_rx_sig(segIndex);
% 
%     % input arg : embedded dimension(m), radius, signal                   
%     realApEnData(frameCount)=ApEn(2, 0.2*std(real(partialData)), real(partialData));
%     imagApEnData(frameCount)=ApEn(2, 0.2*std(imag(partialData)), imag(partialData));
% end
% 
% % load reference data
% load refData_apen.mat;  
%  
% if size(realApEn,1) < limitation
%     limitation = size(realApEn,1) ;
% end
% 
% %% except outlier 
% realApEn = dataConditioning(realApEn);
% imagApEn = dataConditioning(imagApEn);
% 
% %% normalize feature data
% [realApEn, normalizedRealApEnData]= normalization( realApEn, realApEnData );
% [imagApEn, normalizedImagApEnData] = normalization( imagApEn, imagApEnData );  
% 
% %% plot
% scrsz = get(0,'ScreenSize');
% % % [left, bottom, width, height]:
% figure('Name','Approximate entropy','Position',[scrsz(3)/2 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]);
% 
% estimatedSNR=1;
% for floop = 1 : 2
%     ah(floop) = subplot(1,2,floop);  
%     hold(ah(floop),'on');
%     grid(ah(floop),'on');
%     xlabel(ah(floop),'Number of Mode. scheme');
%     ylabel(ah(floop),fSet{floop}); 
%     title((['SNR=',num2str(snr_dB(estimatedSNR)),'dB']));
% end
% 
% boxplot(realApEn(1:limitation,modIdx,estimatedSNR),'parent',ah(1));        
% boxplot(imagApEn(1:limitation,modIdx,estimatedSNR),'parent',ah(2));
% 
% xx = repmat(1:12,size(normalizedRealApEnData,1),1);
% plot(xx,repmat(normalizedRealApEnData(:,estimatedSNR),1,12), 'g.','parent',ah(1));
% xx = repmat(1:12,size(normalizedImagApEnData,1),1);
% plot(xx,repmat(normalizedImagApEnData(:,estimatedSNR),1,12), 'g.','parent',ah(2));
% 
% %% updata result
% set(handles.text_FeaturePoolToFusionCenter,'FontWeight','bold','ForegroundColor',[1 0 0]);
% 
% %% Save data
% 
% handles.realApEn =  realApEn;
% handles.imagApEn =  imagApEn;  
% 
% handles.normalizedRealApEnData = normalizedRealApEnData;
% handles.normalizedImagApEnData = normalizedImagApEnData;
% guidata(hObject,handles);
% 
% set(handles.pushbutton_apen,'BackgroundColor',[222 235 250]/255);


% function pushbutton_instantFreqBasedHist_Callback(hObject, eventdata, handles)
% 
% % parse data from handles
% f_frequency = handles.f_frequency;
% numHist = handles.numHist;
% 
% [binCount, binCenter]=hist(f_frequency, numHist);
% 
% % estimate the number of peaks in histogram
% % threshold = mean(binCount) +  sqrt(var(binCount));
% 
% binCount = binCount ./ max(binCount);
% Param.MinimumHeight=1;
% [NumPeaks, estimatedPeakDistance, threshold]=...
%     peakDetection_140929(binCount.',100,Param);
% %plot graph
% switch NumPeaks
%     case 1
%         set(handles.text_FM,'FontWeight','bold','ForegroundColor',[1 0 0 ]);
%     case {2}
%         set(handles.text_2FSK,'FontWeight','bold','ForegroundColor',[1 0 0 ]);
%     case {3,4}
%         set(handles.text_4FSK,'FontWeight','bold','ForegroundColor',[1 0 0 ]);
%     otherwise        
% end
%  
% set(handles.pushbutton_instantFreqBasedHist,'BackgroundColor',[222 235 250]/255);
% 
% % plot grpah
% cla(handles.axes_monitoring); legend(handles.axes_monitoring,'off')
% % hist(handles.axes_monitoring,f_frequency, numHist);
% bar(handles.axes_monitoring,binCenter, binCount);
% xlim(handles.axes_monitoring,'auto');
% ylim(handles.axes_monitoring,'auto');
% xlabel(handles.axes_monitoring,'Instantaneous Frequency');
% ylabel(handles.axes_monitoring,'Frequency');
% 
% hold(handles.axes_monitoring, 'on');
% plot(binCenter, threshold*ones(1,length(binCenter)),'r');

function pushbutton_makeCHeader_Callback(hObject, eventdata, handles)
%% make C header file
makecheader_150407(handles);

function edit_viewBW_Callback(hObject, eventdata, handles)
function edit_viewBW_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewSNR_Callback(hObject, eventdata, handles)
function edit_viewSNR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewCFO_Callback(hObject, eventdata, handles)
function edit_viewCFO_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewSymRate_Callback(hObject, eventdata, handles)
function edit_viewSymRate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewToneSpacing_Callback(hObject, eventdata, handles)
function edit_viewToneSpacing_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewNumSym_Callback(hObject, eventdata, handles)
function edit_viewNumSym_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_viewSPS_Callback(hObject, eventdata, handles)
function edit_viewSPS_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_ddcDecimationRate_Callback(hObject, eventdata, handles)
%% Parsing
ddcDecimationRate = str2double(get(hObject,'String'));
sam_freq    = handles.sam_freq;
FFTlength   = handles.FFTlength;
% set sampling frequency
sam_freq    = 140e6/ddcDecimationRate;
% set decimation factor
contents  = cellstr(get(handles.popupmenu_decimationFactor,'String'));
num = 1: length(contents);
str = num2str(ddcDecimationRate);
for i = 1 : length(contents)
    if strcmp( contents{i},  str);
        val = i;
    end
end
freqResolution = sam_freq / FFTlength;
%% update GUI states
set(handles.edit_samFreq,'String',num2str(sam_freq));
set(handles.popupmenu_decimationFactor,'Value',val);
set(handles.edit_freqResolution,'String',num2str(freqResolution));
%% save handle structure
handles.previousDecimationFactor = ddcDecimationRate;
handles.sam_freq                 = sam_freq;
guidata(hObject,handles);
function edit_ddcDecimationRate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y = nanmedian(x,dim)
% FORMAT: Y = NANMEDIAN(X,DIM)
% 
%    Median ignoring NaNs
%
%    This function enhances the functionality of NANMEDIAN as distributed
%    in the MATLAB Statistics Toolbox and is meant as a replacement (hence
%    the identical name).  
%
%    NANMEDIAN(X,DIM) calculates the mean along any dimension of the N-D
%    array X ignoring NaNs.  If DIM is omitted NANMEDIAN averages along the
%    first non-singleton dimension of X.
%
%    Similar replacements exist for NANMEAN, NANSTD, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MEDIAN

% -------------------------------------------------------------------------
%    author:      Jan Glscher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.2 $ $Date: 2007/07/30 17:19:19 $

if isempty(x)
	y = [];
	return
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1;
	end
end

siz  = size(x);
n    = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);


% force NaNs to bottom of each column
x = sort(x,1);

% identify and replace NaNs
nans = isnan(x);
x(isnan(x)) = 0;

% new dimension of x
[n m] = size(x);

% number of non-NaN element in each column
s = size(x,1) - sum(nans);
y = zeros(size(s));

% now calculate median for every element in y
% (does anybody know a more eefficient way than with a 'for'-loop?)
for i = 1:length(s)
	if rem(s(i),2) & s(i) > 0
		y(i) = x((s(i)+1)/2,i);
	elseif rem(s(i),2)==0 & s(i) > 0
		y(i) = (x(s(i)/2,i) + x((s(i)/2)+1,i))/2;
	end
end

% Protect against a column of NaNs
i = find(y==0);
y(i) = i + nan;

% permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);

% $Id: nanmedian.m,v 1.2 2007/07/30 17:19:19 glaescher Exp glaescher $



% function pushbutton_spectrumMatching_Callback(hObject, eventdata, handles)
% try
%     updateLoggingMsg('While processing spectrum mathcing , please wait ...',handles);
% 
%     % Parse the data
%     BWSampleIdx         = handles.BWSampleIdx;
%     smoothedPxx_dB      = handles.smoothedPxx_dB;
%     spectrumShapeLimit  = handles.spectrumShapeLimit;
%     
%     rxSpectrumShape  = smoothedPxx_dB(BWSampleIdx);
%     rxSpectrumShape  = rxSpectrumShape(:);
%     
%     %% smoothing
%     rxSpectrumShape = smooth(rxSpectrumShape,3,'moving'); 
% 
%     shapeLength = length(rxSpectrumShape);
%     shapeVector = linspace(1,spectrumShapeLimit,shapeLength);
% 
%     %% interpolate the specturm shape with the length of 256    
%     rxSpectrumShapeInterpolated = interp1(shapeVector,rxSpectrumShape,1:spectrumShapeLimit);
% 
%     %% normalization
%     rxSpectrumShapeNormalized = (rxSpectrumShapeInterpolated- min(rxSpectrumShapeInterpolated)) ./ ...
%         (max(rxSpectrumShapeInterpolated) - min(rxSpectrumShapeInterpolated));
% 
%     %% save data
%     handles.spectrumShape     = rxSpectrumShapeNormalized;
% 
%     %% plot graph  
%     cla(handles.axes_monitoring,'reset');
%     legend(handles.axes_monitoring,'off');
%     
%     plot( rxSpectrumShapeNormalized,'parent',handles.axes_monitoring);
%     xlabel('Samples');  ylabel('Magnitude (dB)'); xlim([1 spectrumShapeLimit]);
%     title('The shape of spectrum','parent',handles.axes_monitoring);
%     %% update GUI states
%     set(handles.pushbutton_spectrumMatching,'BackgroundColor',[222 235 250]/255);
%     set(handles.pushbutton_coarseSymRateEstimation,'Enable','on');
% 
%     guidata(hObject,handles);
%     updateLoggingMsg('Complete spectrum shape matching section ',handles);
% catch err
%     rethrow(err);
%     updateLoggingMsg('Error in spectrum shape matching section ',handles);
% end


% --- Executes during object deletion, before destroying properties.
function fig_AMRtestbed_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to fig_AMRtestbed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
