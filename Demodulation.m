function varargout = Demodulation(varargin)
% DEMODULATION MATLAB code for Demodulation.fig
%      DEMODULATION, by itself, creates a new DEMODULATION or raises the existing
%      singleton*.
%
%      H = DEMODULATION returns the handle to a new DEMODULATION or the handle to
%      the existing singleton*.
%
%      DEMODULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMODULATION.M with the given input arguments.
%
%      DEMODULATION('Property','Value',...) creates a new DEMODULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Demodulation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Demodulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%  pushbutton_viewWaveform_Callback(hObject, eventdata, handles)

% Edit the above text to modify the fresponse to help Demodulation

% Last Modified by GUIDE v2.5 12-May-2015 18:59:17
% Last Modified by hand  15.03.03

% Begin initialization code - DO NOT EDIT

% Reference
%
% [1] Mengali, Umberto, and Aldo N. D'Andrea, SynchronizationTechniques for 
%     Digital Receivers, New York, Plenum Press,1997. 
% [2] Michael Rice, "Digital Communications - A Discrete-Time Approach".
%
% Functions
% pushbutton_exeDemod_Callback
% function demodSig = PSKdemodulation(handles)
% function demodSig = FSKdemodulation(handles)
% 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Demodulation_OpeningFcn, ...
                   'gui_OutputFcn',  @Demodulation_OutputFcn, ...
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

function Demodulation_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for Demodulation
handles.output          = hObject;
handles.initialPhase    = 0;
handles.carrierOffset   = 0;
handles.decoder         = 'CVSD';
handles.filterType      = 'Raised Cosine';
set(handles.edit_initialPhase,'String',num2str(handles.initialPhase));

%% Update handles structure
guidata(hObject, handles);

function varargout = Demodulation_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function pushbutton_input_Callback(hObject, eventdata, handles)
       
[filename,filepath]=uigetfile({'*.*','All Files'},...
    'Select Data File 1');
if filename ~= 0

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
            sigLenLimitEnd      = floor(length( initial_uint32_Data ));       
            initial_uint32_Data = initial_uint32_Data(sigLenLimitStart:sigLenLimitEnd);

            % Discard first 16byte
            uint32_Data = initial_uint32_Data(5:end);
            clear initial_uint32_Data;

            sigLen = length(uint32_Data ); 
            % Discard remaining byte
            uint32_Data = uint32_Data(1:(sigLen-mod(sigLen,4))); 
            % Flip data by 16byte 
            temp = reshape(uint32_Data,4,[]);
            uint32_Data = temp(4:-1:1,:);
            clear temp;
%                 uint32_Data = temp2;           

            % Find TOA Indexes
            idxTOA = find( bitand(uint32_Data,hex2dec('c0000000')) == hex2dec('c0000000'));

            % discard samples located in before the first TOA
            % uint32_Data=initial_uint32_Data(idxTOA(1): end);
            % uint32_Data = initial_uint32_Data;
%                 idxSet = true(1,length(uint32_Data));

            % Discard 5 of samples in the next of TOA
            deleteIdx= [];
            for ioop = 1 : length(idxTOA)
            %    deleteIdx = [deleteIdx idxTOA(ioop) : idxTOA(ioop)+ 5]; 
                 deleteIdx = [deleteIdx idxTOA(ioop) : (idxTOA(ioop)+ 5) ]; % edited in 150127
            end
            uint32_Data( unique(deleteIdx) ) = [];
    %%%%%%%%%%%%%%%%%%%%%%%% IQ data structure %%%%%%%%%%%%%%%%%%%%%%%%
    % Type      |  header  |        Inphase      |       Quadrature   |     
    % Bit index |  31   30 |     29     ~    15  |   14      ~      0 |           
    %           | TOA : 10 | (sign bit)          | (sign bit)         |
    %           | IQ  : 01 |                     |                    |
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oasis = bitand( uint32_Data,hex2dec('3fff8000'));
            i_tmp = bitshift(oasis,-15); clear oasis;

            q_tmp = bitand( uint32_Data,hex2dec('00007fff'));
            clear uint32_Data;

            % two's complement
            maskingHex4000 = hex2dec('4000'); % sign bit
            % find negative value
            comp_I_Idx = bitand(i_tmp,maskingHex4000) == maskingHex4000;

            IData = zeros(1,length(comp_I_Idx));
            masking = hex2dec('3fff'); 
            % Take two's complement
            tmp = (masking - bitand( i_tmp(comp_I_Idx),masking )) +1; 
            % Save value as double-precision floating-point format in MATLAB
            IData(comp_I_Idx) = -double(tmp).';  
            tmp = bitand(i_tmp(~comp_I_Idx),masking);
            IData(~comp_I_Idx) = double(tmp);

            clear comp_I_Idx;
            comp_Q_Idx = bitand(q_tmp,maskingHex4000) == maskingHex4000;
            QData = zeros(1,length(comp_Q_Idx));
            tmp = (masking - bitand( q_tmp(comp_Q_Idx),masking )) +1;
            QData(comp_Q_Idx) =  -double(tmp).';
            tmp = bitand(q_tmp(~comp_Q_Idx),masking);
            QData(~comp_Q_Idx) = double(tmp);
            clear comp_Q_Idx; clear q_tmp;

            sigLenLimitEnd = length(IData);
            sigLenLimitStart = 1;
            IQData = IData + 1i*QData;    % IQ data 
            clear IData; clear QData;
            % remove mean of IQData
            % IQData = IQData-mean(IQData);
            % make average power to 1
            % IQData = sqrt(length(IQData)) * IQData / norm(IQData,2);                 
            rx_sig = IQData;    
        else
            updateLoggingMsg('Move the exe file to the folder including diq data',handles);

        end
    case '.mat'
        load(filename);
        
        rx_sig                           = AMRtestbedHandles.rx_sig;
        handles.modType                  = AMRtestbedHandles.coarseModType;
        handles.sam_freq                 = AMRtestbedHandles.sam_freq;                   % sampling frequency
        handles.bandwidth                = AMRtestbedHandles.coarseBandWidth;
        handles.carrierOffset            = AMRtestbedHandles.coarseFreqOffset;
        handles.samplesPerSymbol         = AMRtestbedHandles.coarseSamplesPerSymbol;
        handles.symbolRate               = AMRtestbedHandles.coarseSymbolRate;
        handles.toneSpacing              = AMRtestbedHandles.coarseToneSpacing;
        handles.modulationOrder          = AMRtestbedHandles.modulationOrder;
        handles.decimationFactor         = AMRtestbedHandles.previousDecimationFactor;
        handles.coarseModType            = AMRtestbedHandles.coarseModType;
        handles.defaultRx_sig            = rx_sig;
              
        switch handles.modType 
            case 'AM',                            modType = 1;
            case 'FM',                            modType = 2;
            case {'ASK','2ASK'},                  modType = 3;
            case {'FSK','2FSK','4FSK','8FSK'},    modType = 4;       
            case {'PSK','BPSK', 'QPSK', '8PSK'},  modType = 5;       
            case {'QAM','16QAM','32QAM','64QAM'}, modType = 6;       
        end

        sigLenLimitEnd = length(rx_sig);
        sigLenLimitStart = 1;     
        % AM,FM,2ASK,FSK,PSK,QAM
        set(handles.popupmenu_modType,'Value',modType);
        set(handles.edit_viewStartIdx,'String',num2str( sigLenLimitStart));    
        set(handles.edit_viewEndIdx,'String',num2str( sigLenLimitEnd));    
        set(handles.edit_numberOfInputSamples,'String',num2str(sigLenLimitEnd));
        set(handles.edit_modulationOrder,'String',num2str(handles.modulationOrder));
        set(handles.edit_samFreq,'String',num2str(handles.sam_freq));
        set(handles.edit_carrierOffset,'String',num2str(handles.carrierOffset));
        set(handles.edit_symbolRate,'String',num2str(handles.symbolRate));
        set(handles.edit_toneSpacing,'String',num2str(handles.toneSpacing));                
        set(handles.popupmenu_resampleRate,'Value',3);    
        set(handles.text_dispFileName,'String','AMRtestbedHandles.mat');

    case '.fig'
        % Open FIG-file with guide command.
        guide(filename)
        rx_sig = [1 2];

    case {'.wav', '.mp3'}
         % When size of a file is too large, 
         % restrict the number of loaded samples
         [rx_sig, handles.sam_freq] = ...
              audioread(filename,samples);             

        sigLenLimit = length(rx_sig);
        sigLenLimitStart = 1;
        sigLenLimitEnd = sigLenLimit;   

    otherwise
        try
            % Use open for other file types.
            open(filename)
        catch ex
            errordlg(...
              ex.getReport('basic'),'File Type Error','modal')
        end

end

%% update GUI states
set(handles.text_dispFileName,'String',filename);
set(handles.edit_numberOfInputSamples,'String',num2str(sigLenLimitEnd-sigLenLimitStart)); 
set(handles.edit_viewEndIdx,'String',num2str(sigLenLimitEnd));
set(handles.edit_viewStartIdx,'String',num2str( sigLenLimitStart));
set(handles.edit_numberOfInputSamples,'String',num2str(length(rx_sig)));
set(handles.popupmenu_resampleRate,'Value',3);    

%% save data
handles.soundFs       = 44100; 
handles.filename      = filename;
handles.rx_sig        = rx_sig;
handles.defaultRx_sig = rx_sig;
handles.viewStartIdx  = sigLenLimitStart;
handles.viewEndIdx    = sigLenLimitEnd;
guidata(hObject, handles);

%% Update GUI states

end

guidata(hObject,handles);

function edit_samFreq_Callback(hObject, eventdata, handles)
handles.sam_freq = str2double ( get(hObject,'String') );
updateLoggingMsg(['Set sampling frequency :',num2str(handles.sam_freq),' [Hz]'],handles);
guidata(hObject,handles);

function edit_samFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_carrierOffset_Callback(hObject, eventdata, handles)
handles.coarseFreqOffset = str2double ( get(hObject,'String') );
updateLoggingMsg(['Set frequency offset :',num2str(handles.coarseFreqOffset),' [Hz]'],handles);
guidata(hObject,handles);

function edit_carrierOffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_symbolRate_Callback(hObject, eventdata, handles)
handles.symbolRate = str2double ( get(hObject,'String') );
updateLoggingMsg(['Set symbol rate :',num2str(handles.symbolRate),' [Hz]'],handles);
guidata(hObject,handles);

function edit_symbolRate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_toneSpacing_Callback(hObject, eventdata, handles)
handles.toneSpacing = str2double ( get(hObject,'String') );
updateLoggingMsg(['Set tone spacing :',num2str(handles.toneSpacing),' [Hz]'],handles);
guidata(hObject,handles);

function edit_toneSpacing_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_modulationOrder_Callback(hObject, eventdata, handles)
handles.modulationOrder =str2double(get(hObject,'String'));
updateLoggingMsg(['Set modulation order to ',num2str(handles.modulationOrder)],handles);
guidata(hObject,handles);

function edit_modulationOrder_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_initialPhase_Callback(hObject, eventdata, handles)
handles.initialPhase = eval(get(hObject,'String'));
updateLoggingMsg(['Set initial phase to ',num2str(handles.initialPhase),' [rad]'],handles);
guidata(hObject,handles);

function edit_initialPhase_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_filterType_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String')) ;
handles.filterType = contents{get(hObject,'Value')} ;
updateLoggingMsg(['Set filter type :',handles.filterType],handles);
guidata(hObject,handles);

function popupmenu_filterType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_decoder_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String')) ;
handles.decoder = contents{get(hObject,'Value')} ;
updateLoggingMsg(['Set decoder :',handles.decoder],handles);
guidata(hObject,handles);

function popupmenu_decoder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_exeDemod_Callback(hObject, eventdata, handles)

updateLoggingMsg('While demodulation processing ...',handles);
sam_freq                  = handles.sam_freq;
modType                   = handles.modType;
modulationOrder           = handles.modulationOrder;
symbolRate                = handles.symbolRate;
viewStartIdx              = handles.viewStartIdx;
viewEndIdx                = handles.viewEndIdx;
carrierOffset             = handles.carrierOffset;
toneSpacing               = handles.toneSpacing;          
initialPhase              = handles.initialPhase;

samplesPerSymbol          = sam_freq / symbolRate;
% samplesPerSymbol          = round(samplesPerSymbol);
symbolRate                = ceil(symbolRate);
sam_freq                  = ceil(sam_freq);       

sam_period = 1/ sam_freq;
if isfield(handles,'resampledRxSig')
    rx_sig              = handles.resampledRxSig;
else
    rx_sig              = handles.rx_sig;
end

rx_sig = rx_sig(viewStartIdx : viewEndIdx);
%% exploit the instantaneous information 
% envelope=abs(rx_sig);                   % find the envelope 
% myAmp = envelope-mean(envelope);
% 
% nn = 1000;
% xx = (0: nn-1); %/sam_freq;
% xx2 = xx(1:end-1);
% 
% myAmp = myAmp(1:nn);
% myPhase = unwrap(angle(rx_sig(1:nn)));      
% myFreq = diff(myPhase); 
% myFreq = smooth(myFreq,20,'moving');  
% 
% figure('name',handles.filename);
% subplot(4,1,1), plot(xx,real(rx_sig(1:nn)),'-r.',xx,imag(rx_sig(1:nn)));title('IQ plot'); grid on;
% subplot(4,1,2), plot(xx,myAmp,'-b');title('Instantaneous amplitude');grid on;
% subplot(4,1,3), plot(xx,myPhase/10,'-b');title('Instantaneous phase');grid on;
% subplot(4,1,4), plot(xx2,(myFreq),'-r.'); title('Instantaneous frequency');grid on;
% 
%  myPhase = unwrap(angle(rx_sig));
% myFreq = diff(myPhase); 
% myFreq = smooth(myFreq,20,'moving');  
% 
% figure('name',handles.filename);
% subplot(3,1,1), hist(myAmp, 100); title('Histogram of IA')
% subplot(3,1,2), hist(myPhase,100); title('Histogram of IP')
% subplot(3,1,3), hist(myFreq, 100); title('Histogram of IF')
% mean(myFreq)
% a=4;
%% root-raised cosine filtering
% % Shape of the pulse shaping filter
% shape = 'Square Root Raised Cosine';
% Nsym = 8; beta = 0.35;
% % Specifications of the raised cosine filter with given order in symbols
% sqrtRcosSpec = fdesign.pulseshaping(samplesPerSymbol, shape,...
%     'Nsym,beta', Nsym, beta);
% % Design and normalize filter.
% sqrtRcosFlt = design(sqrtRcosSpec);
% normFact = max(sqrtRcosFlt.Numerator);
% sqrtRcosFlt.Numerator = sqrtRcosFlt.Numerator / normFact;
% % Upsample and filter.
% rx_sig = filter(sqrtRcosFlt, rx_sig);

%% demodulation
switch modType
    case 'AM'
        % -----------------------------------------------------------------
        % noncoherent demodulation : envelope detection
        % We use Hilbert transform to find the envelope of the bandpass signal
        % The Hilbert transform function in Matlab, denoted by hilbert.m,
        % generates the analytic signal z(t).
        % The real part of z(t) is the original sequence, and its imaginary part
        % is the Hilbert transform of the original sequence.
        % a : modulation index (0, 1] -->   (1+a*m(t))*cos(wc*t)
        % m_max : maximum amplitude of information signal.
        % -----------------------------------------------------------------
        % complex envelope signal having little frequency offset
        envelope=abs(rx_sig);                   % find the envelope 
%          envelope = envelope / 0.9; %max(envelope);
        demodSig = envelope-mean(envelope);
        % dem1=m_max*(envelope-1)/a; % remove dc and rescale

    case 'FM'
        freqdev = 20E3;  % Default frequency deviation
  
        demodSig = (1/(2*pi*freqdev))*[zeros(1,size(rx_sig,2)); diff(unwrap(angle(rx_sig)))*sam_freq]; % complete C conversion 
        
    case{'ASK', '2ASK'}
        envelope=abs(rx_sig);                   % find the envelope 
        myAmp = envelope-mean(envelope);

       %% Symbol timing recovery 2 - ZCTED, 2 samples per symbol
        %Down-sampling
        D = floor(samplesPerSymbol / 2);
        x = myAmp(1:D:end);
        
       %% Param
        N = 2;  % samples per symbol
        noiseBandwidth = 0.01; % BnTs, normalized to the bit rate
        dampingFactor = 1; %1/sqrt(2); % Zeta
        phaseDetectorGain = 0.6; % kp % time to approach the steady state
        k0 = -1; % the controller is a decrementing modulo-1 counter
        const = noiseBandwidth/(dampingFactor + (1/(4*dampingFactor)));
        
        den = (const*4*dampingFactor/N);
        nom = (1+ const*2*dampingFactor/N + (const/N)^2);
        
        k1 = (den/nom) /(phaseDetectorGain*k0);
    
        den2 = (const*const*4/(N*N));
 
        k2 = (den2/nom) /(phaseDetectorGain*k0);
        TEDBuff = 0;    CNT_next = 1;   mu_next = 0; %initial timing offset
        underflow = 0;  old_underflow =1;
          k = 1;  xx=[];
              vi = 0;     % initial state of loop filter  
        for n=4:length(x) % for causal system
            %evaluate arithmetic expressions in topological order
            CNT = CNT_next; % etha
            mu = mu_next;
            v2 = 1/2*[1, -1, -1, 1]*x(n:-1:n-3);  % Farrow structure for the piecewise parabolic
            v1 = 1/2*[-1, 3, -1, -1]*x(n:-1:n-3); % interpolator, 알파로 묶고 정리하면 현재식이 된다
            v0 = x(n-2); % baseindex
            xI = (mu*v2 + v1)*mu + v0;                
            if underflow == 1
                xx(k) = xI;
                k = k + 1;               
            end

            if underflow == 1 && old_underflow == 0
                
                rexI = real(xI);                imxI = imag(xI);
                reTEDBuff = real(TEDBuff);      imTEDBuff = imag(TEDBuff);
                % Zero-crossing TED
                e = reTEDBuff(1) * (sign(reTEDBuff(2)) - sign(rexI)) + ...
                    imTEDBuff(1) * (sign(imTEDBuff(2)) - sign(imxI));
            else
                e = 0;
            end
            vp = k1*e;          % proportional component of loop filter
            vi = vi + k2*e;     % integrator componenet of loop filter
            v = vp + vi;        % loop filter output
            W = 1/N + v;        % NCO control word

            % update registers
            if underflow == 0 && old_underflow == 0
                TEDBuff = TEDBuff;                          % skip current sample
            elseif underflow == 0 && old_underflow == 1
                TEDBuff = [xI; TEDBuff(1)];                 % normal operation
            elseif underflow ==1 && old_underflow == 0
                TEDBuff = [xI; TEDBuff(1)];                 % normal operation
            elseif underflow ==1 && old_underflow == 1
                TEDBuff = [xI; 0];                          % stuff missing sample
            end
            
            CNT_next = CNT - W;             % update counter value for next cycle
            if CNT_next <0                  % test to see if underflow has occurred
                CNT_next = 1 + CNT_next;    % reduce counter value modulo-1 if underflow
                old_underflow = underflow;  
                underflow = 1;              % set underflow flag
                mu_next = CNT/W;            % update mu
            else
                old_underflow = underflow;
                underflow = 0;
                mu_next = mu;
            end    
        end
        myAmp = xx;
        scatterplot(xx);
        
        demodSig = zeros(1,length(myAmp));
        highIdx = myAmp > 0 ;
        demodSig(highIdx) = 1;
        demodSig(~highIdx) = 0;

%         highIdx = myAmp(1:samplesPerSymbol:end) > 0 ;
%         demodSig = zeros(1,length(myAmp(1:samplesPerSymbol:end)));
%         demodSig(highIdx) = 1;
%         demodSig(~highIdx) = 0;
        
    case {'FSK','2FSK','4FSK','8FSK'}
        demodSig = FSKdemodulation(handles);
    case {'PSK','BPSK','QPSK','8PSK'} 
        demodSig = PSKdemodulation(handles);
    case {'QAM','32QAM','16QAM','64QAM'}
        % De-rotate
        rx_sig = rx_sig .* exp(-1i*ini_phase);

        % Precompute for later use
        sqrtM = sqrt(modulationOrder);

        % Inphase/real rail
        % Move the real part of input signal; scale appropriately and round the
        % values to get index ideal constellation points
        rIdx = round( ((real(rx_sig) + (sqrtM-1)) ./ 2) );
        % clip values that are outside the valid range 
        rIdx(rIdx <= -1) = 0;
        rIdx(rIdx > (sqrtM-1)) = sqrtM-1;

        % Quadrature/imaginary rail
        % Move the imaginary part of input signal; scale appropriately and round 
        % the values to get index of ideal constellation points
        iIdx = round(((imag(rx_sig) + (sqrtM-1)) ./ 2));
        % clip values that are outside the valid range 
        iIdx(iIdx <= -1) = 0;
        iIdx(iIdx > (sqrtM-1)) = sqrtM-1;

        % compute output from indices of ideal constellation points 
        demodSig = sqrtM-iIdx-1 +  sqrtM*rIdx;
    otherwise
end
        %%-------------------------------------------
        % Noncoherent FSK Demodulation, 
        % -------------------------------------------
        % refer to ex10_3.m
        % Theoretical envelope detection        
%         dr_dt = diff(xm);
%         z = hilbert(dr_dt);         % get analytic signal 
%         envelope = abs(z);            % find the envelope 
%         demod = envelope;
%         subplot(211), waveform(x, fs);
%         title('information signal');
%         subplot(212), cw_waveform(demod,fs); 
%         title('envelope detector output');
%         r = axis; range = [r(1) r(2) 0 1];
%         axis(range);
updateLoggingMsg('Demodulation complete',handles);

%% update handles structure
handles.viewStartIdx = 1;
handles.viewEndIdx   = length(demodSig);
handles.sam_freq     = sam_freq;
handles.rx_sig       = demodSig;

%% temporary coding 150317
% handles.rx_sig = rx_sig;

guidata(hObject,handles);

function demodSig = FSKdemodulation(handles)

try
    sam_freq          = handles.sam_freq;
    modulationOrder   = handles.modulationOrder;
    symbolRate        = handles.symbolRate;
    viewStartIdx      = handles.viewStartIdx;
    viewEndIdx        = handles.viewEndIdx;
    sam_freq          = ceil(sam_freq);       
    toneSpacing       = handles.toneSpacing;
    
    if isfield(handles,'resampledRxSig')
        rx_sig              = handles.resampledRxSig;
    else
        rx_sig              = handles.rx_sig;
    end

    rx_sig = rx_sig(viewStartIdx : viewEndIdx);
    phaseOffset = 0;
    samplesPerSymbol = floor(sam_freq/symbolRate);
    %%-------------------------------------------
    % Noncoherent FSK Demodulation
    % -------------------------------------------       
    % please refer to fskdemod.m
    % Check that the maximum transmitted frequency does not exceed Fs/2
    maxFreq = ((modulationOrder-1)/2) * toneSpacing;
    if (maxFreq > sam_freq/2)
        error(message('comm:fskdemod:maxFreq'));
    end
    % Preallocate memory
    demodSig = zeros(1, floor(length(rx_sig)/samplesPerSymbol));
    %         nSamp = length(demodSig);
    % Define the frequencies used for the demodulator.  
    freqs = (-(modulationOrder-1)/2 : (modulationOrder-1)/2) * toneSpacing;
    % Use the frequencies to generate M complex tones which will be multiplied with
    % each received FSK symbol.  The tones run down the columnns of the "tones"
    % matrix.
    t = [0 : 1/sam_freq : samplesPerSymbol/sam_freq - 1/sam_freq]';        
    phase = 2*pi*t*freqs + 2*pi*phaseOffset;
    tones = exp(-1i*phase);            

    for iSym = 1 : length(rx_sig)/samplesPerSymbol

        % Load the samples for the current symbol
        yTemp = rx_sig( (iSym-1)*samplesPerSymbol+1 : iSym*samplesPerSymbol);

        % Replicate the received FSK signal to multiply with the M tones
        yTemp = yTemp(:, ones(modulationOrder,1));

        % Multiply against the M tones
        yTemp = yTemp .* tones;
    %             figure; plot(real(yTemp(:,1)),'-b.'); hold on; plot(imag(yTemp(:,1)),'-r.');
    %             figure; plot(real(yTemp(:,2)),'-b.'); hold on; plot(imag(yTemp(:,2)),'-r.');

        % Perform the integrate and dump, then get the magnitude.  Use a
        % subfunction for the integrate and dump, to omit the error checking.
        yMag = abs(intanddump(yTemp, samplesPerSymbol));
    %             yTemp = mean(reshape(yTemp, Nsamp, length(rx_sig)/samplesPerSymbol), 1);           
    %             yMag = abs(yTemp);

        % Choose the maximum and assign an integer value to it.  Subtract 1 from the
        % output of MAX because the integer outputs are zero-based, not one-based.
        [~, maxIdx] = max(yMag, [], 2);

        demodSig(iSym) = maxIdx-1;
    end

    switch modulationOrder
        case 4
            mappingTable =[0 0 1 1;...
                           0 1 0 1];
            bitStream = mappingTable(:,demodSig);
            demodSig = bitStream(:);
        case 8
            mappingTable =[0 0 0 0 1 1 1 1 ; ...
                           0 0 1 1 0 0 1 1 ; ...
                           0 1 0 1 0 1 0 1];
            bitStream = mappingTable(:,demodSig);
            demodSig = bitStream(:);       
    end

catch err
    rethrow(err);
end

function demodSig = PSKdemodulation(handles)
%%
try
sam_freq          = handles.sam_freq;
modulationOrder   = handles.modulationOrder;
symbolRate        = handles.symbolRate;
viewStartIdx      = handles.viewStartIdx;
viewEndIdx        = handles.viewEndIdx;
sam_freq          = ceil(sam_freq);       
sam_period        = 1/ sam_freq;
initialPhase      = handles.initialPhase;
isDQPSK = 1;

if isfield(handles,'resampledRxSig')
    rx_sig              = handles.resampledRxSig;
else
    rx_sig              = handles.rx_sig;
end

rx_sig = rx_sig(viewStartIdx : viewEndIdx);
%% Automatic Gain Control
x = rx_sig;

updateLoggingMsg('While AGC processing ...',handles);
%% AGC function
% Default parameter
UpdatePeriod = 100;
StepSize = 0.1;
MaximumGain = 30; 
ReferenceLevel = 1;
Gain = 1;

% re-assign
g = Gain;    
agcOut = zeros(size(x));
K = StepSize;    
ref = ReferenceLevel;   
maxGain = MaximumGain;
minGain = eps; 
updatePeriod = UpdatePeriod;
numSubFrames = length(x)/updatePeriod;

% Linear method
for p=0:numSubFrames-1
    indices = p*updatePeriod + (1:updatePeriod);
    agcOut(indices) = x(indices) * g;
    z = mean(abs(agcOut(indices)));  % same as  z = rectifier(y(indices));
    e = ref - z;
    g = g + K*e;
    if g < minGain
        g = minGain;
    elseif g > maxGain
        g = maxGain;
    end
end

%% Coarse frequency offset compensation
updateLoggingMsg('While coarse frequency offset estimation processing ...',handles);

order = 4;
nonLinearity = agcOut.^order;
desireFreqResolution = (symbolRate/100) / 10;
NFFT = 512;
for i=1:500
   NFFT= NFFT * 2;
   if (sam_freq/desireFreqResolution) < NFFT
        break; 
   end
end

spectrum = ( fft(nonLinearity,NFFT) );
absSpectrum = abs(spectrum);
% freq = (0 : NFFT -1)/NFFT * sam_freq; 
freq = (-NFFT/2 : NFFT/2 -1)/NFFT * sam_freq; 
% figure; plot(freq,fftshift(absSpectrum));
% 
%         searchRange = freq <= 2e3; % set the searching range
[~, indexFreqOffset ] = max(absSpectrum);

% De-rotation
if isDQPSK
    %     http://cp.literature.agilent.com/litweb/pdf/ads2008/timed/ads2008/DQPSK_Pi4DemodSync.html
%     figure; plot(fftshift(absSpectrum));
    symRateLineSpec = freq(indexFreqOffset);
    estimatedFreqOffset = -1170 / 4; %(symRateLineSpec - symbolRate/2)/4;
    tt= (0:(length(agcOut)-1)) /sam_freq; tt= tt(:);
    coarseFreqCompSignal = agcOut.*exp(-1i*2*pi*estimatedFreqOffset*tt);

%     absSpectrum = abs(fftshift( fft(coarseFreqCompSignal.^4,NFFT) ));
%     figure; plot(freq,fftshift(10*log10(absSpectrum)));
else
    estimatedFreqOffset = freq(indexFreqOffset)/order;
    tt= (0:(length(agcOut)-1)) /sam_freq; tt= tt(:);
    coarseFreqCompSignal = agcOut.*exp(-1i*2*pi*estimatedFreqOffset*tt);
    
end

%% Sample-rate conversion by a factor of P/Q (used to
% downconvert from xxkHz to xxkHz)
% See Chaning singal Sampling Rate example
updateLoggingMsg('While sample-rate conversion processing ...',handles);

symbolRate =round(symbolRate);
desireFs = symbolRate * 2;
%  
[P,Q] = rat( desireFs / sam_freq);   % Detemine the Interpolation/decimation factors.
abs(P/Q*sam_freq-desireFs)
% 
% % L  = 147; M = 160;                 
Hm = mfilt.firsrc(P,Q);              % We use the default filter
sampleConversionSig = filter(Hm,coarseFreqCompSignal);                    %   

% dFactor = floor( (sam_freq/symbolRate)/2 );
% sampleConversionSig = downsample(coarseFreqCompSignal,dFactor);
updateLoggingMsg('While fine carrier frequency offset and symbol-timing recovery processing ...',handles);
%% Paramter for fine carrier frequency offset  - in [2] p105
K = 1;
A = 1/sqrt(2);
% Look into model for details for details of PLL parameter choice. 
% Refer equation 7.30 in [2] 
% K_p for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM),
% QPSK could be treated as two individual binary PAM
% PhaseErrorDetectorGain = 2*K*A^2+2*K*A^2; % QPSK
PhaseErrorDetectorGain = 1; % 8PSK
% PhaseRecoveryGain = 1; % K_0 for Fine Frequency Compensation PLL
PhaseRecoveryGain = 2; % equal to the number of samples per symbol

PhaseRecoveryLoopBandwidth = 0.005; % Normalized loop bandwidth for fine frequency compensation
PhaseRecoveryDampingFactor = 1; % Damping Factor for fine frequency compensation
PostFilterOversampling = 2;
% Refer C.57 to C.61 in [2] for K1 and K2
theta = PhaseRecoveryLoopBandwidth/...
    (PhaseRecoveryDampingFactor + ...
    0.25/PhaseRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*PhaseRecoveryDampingFactor*theta + theta*theta;
K1 = (4*PhaseRecoveryDampingFactor*theta/d)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);
K2 = (4*theta*theta/d)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);

IntegratorGain = K2;
ProportionalGain = K1;

DigitalSynthesizerGain = -1;
pPhase = 0;

%% Parameter-Timing recovery
% Parameters
N = 2;  % samples per symbol
noiseBandwidth = 0.01; % BnTs, normalized to the bit rate
dampingFactor = 1; %1/sqrt(2); % Zeta
phaseDetectorGain = 2.7; % kp % time to approach the steady state
k0 = -1; % the controller is a decrementing modulo-1 counter
const = noiseBandwidth/(dampingFactor + (1/(4*dampingFactor)));
den = (const*4*dampingFactor/N);
nom = (1+ const*2*dampingFactor/N + (const/N)^2);
ProportionalGainTiming = (den/nom) /(phaseDetectorGain*k0);

den2 = (const*const*4/(N*N));
IntegratorGainTiming = (den2/nom) /(phaseDetectorGain*k0);

pAlpha = 0.5;
pDelayBuffer1 = 0; pDelayBuffer2 = 0; pDelayBuffer3 = 0;    
pTEDDelay1 = 0; pTEDDelay2 = 0; 
pDelayStrobe = 0;

pNCODelay = 0;
pRegTemp = 0;
pMU = 0;
pStrobe = 0;
timingLoopOut = 0; pCount = 1; % timing loop output 
fineCompSignalBuff = 0;
loopFiltOut=0; DDSOut=0; loopFiltOutTiming = 0.0;
figure; axis([ -1.5 1.5 -1.5 1.5]);
for i=1:length(sampleConversionSig)
    %% Carrier offset recovery  - 2 samples per symbol   
    % Complex phase shift
    fineCompSignal = sampleConversionSig(i) * exp(1j*pPhase);
    fineCompSignalBuff(i) = fineCompSignal;
    
    % Find phase error
    re = real(fineCompSignal); im = imag(fineCompSignal);
%     plot(re,im,'r.'); hold on; plot(real(sampleConversionSig(i)),imag(sampleConversionSig(i)),'b.');
%     pause(0.0001);
    if isDQPSK        
        %for 8PSK
        % http://kr.mathworks.com/help/comm/ref/comm.carriersynchronizer-class.html#bulwrmb
        a = sqrt(2)-1;
        if re >= im
            phErr =  a*sign(re) *im - sign(im)*re; 
        else
            phErr = sign(re)*im - a*sign(im)*re;
        end
    else
        % QAM and QPSK
        phErr = sign(re)*im - sign(im)*re;
    end
    
    loopFiltOut = loopFiltOut +(phErr*IntegratorGain);        % loop filter output

    % Direct Digital Synthesizer
    DDSOut = DDSOut + (phErr*ProportionalGain + loopFiltOut); % pIntegrator
    pPhase =  DigitalSynthesizerGain * DDSOut;

    %% Timing recovery  - ZCTED, 2 samples per symbols
    % Parabolic interpolator - Farrow structure, p 471 in [1]
    K = -pAlpha;      
    interpFiltOut = pDelayBuffer2 + ...       
        pMU*(K*(fineCompSignal+pDelayBuffer2+pDelayBuffer3)+...
        (1-K)*pDelayBuffer1)...
        + K*(pDelayBuffer1+pDelayBuffer2-fineCompSignal-pDelayBuffer3)...
        *pMU^2;

    % Update delay buffers
    pDelayBuffer3 = pDelayBuffer2;
    pDelayBuffer2 = pDelayBuffer1;
    pDelayBuffer1 = fineCompSignal;    

    % Timing Error Detector
    % e = TED(obj,interpFiltOut);
   if pStrobe && pDelayStrobe~=pStrobe
        e = real(pTEDDelay1) * (sign(real(pTEDDelay2)) - ...
            sign(real(interpFiltOut))) + imag(pTEDDelay1) * ...
            (sign(imag(pTEDDelay2)) - sign(imag(interpFiltOut)));
    else
        e = 0;
    end     

    if pDelayStrobe~=pStrobe
        % Shift contents in delay register
        pTEDDelay2 = pTEDDelay1;
        pTEDDelay1 = interpFiltOut;
    elseif pStrobe
        % Two consecutive high strobes
        pTEDDelay2 = 0; % Stuff missing sample
        pTEDDelay1 = interpFiltOut;
    end
    % If both current and previous enable signals are 0, skip current sample
    % and keep the delayed signals unchanged. (Bit stripping)
    pDelayStrobe = pStrobe;   

    % Loop filter
    loopFiltOutTiming = loopFiltOutTiming + IntegratorGainTiming*e;     % integrator componenet of loop filter

    % Updates timing difference for Interpolation Filter
    % [underflow,pMU] = NCO_control(obj,e,loopFiltOut);...
    Delta = e*ProportionalGainTiming + loopFiltOutTiming; % If loop is in lock, Delta would be small
    modVal = mod(pNCODelay,1);
    cc = Delta+1/PostFilterOversampling;
    counter = modVal-cc; % decrementing counter
    if counter < 0
        underflow = 1;
        pRegTemp = modVal/cc;
    else
        underflow = 0;
    end 

    pMU = pRegTemp;
      % Update delay buffer
    pNCODelay = counter; %mod(pNCODelay,1)-Delta-1/PostFilterOversampling;  
    if pStrobe>0
       timingLoopOut(pCount) = interpFiltOut;
%          plot(real(timingLoopOut(pCount)),imag(timingLoopOut(pCount)),'g.');
       pCount = pCount + 1;
      
    end

    pStrobe = underflow;   

end
rx_sig = timingLoopOut;
scatterplot(fineCompSignalBuff(1:4000));
title('After fine freq. offset estimation processing');

scatterplot(rx_sig);
title('After timing error estimation processing');
%% symbol timing recovery
%         samplesPerSymbol = floor(samplesPerSymbol);
%         
% %         hSync = comm.GardnerTimingSynchronizer('SamplesPerSymbol', samplesPerSymbol, ...
% %         'ErrorUpdateGain', 0.07);
%     
%         hSync = comm.EarlyLateGateTimingSynchronizer(...
%                           'SamplesPerSymbol',samplesPerSymbol, ...
%                           'ErrorUpdateGain', 0.075);
%         limit = floor( length(rx_sig) / samplesPerSymbol );
%         rx_sig = rx_sig(1:limit*samplesPerSymbol);
%         
% %         figure; plot(real(rx_sig),'-r.'); hold on; plot(imag(rx_sig),'-b.');
%         [rx_sig, phase] = step(hSync, rx_sig);
%         figure; plot(real(rx_sig2),'-k.'); hold on; plot(imag(rx_sig2),'-g.');
%         scatterplot(rx_sig);
%         eyediagram(rx_sig(1:1000),samplesPerSymbol*2)
% % 
        %% method 3 - Approximate ML Estimation
%         order = 2;
%         NFFT = 2^(nextpow2(length(rx_sig))-1);
%         ff = (-NFFT/2 : NFFT/2-1)/ NFFT * sam_freq;
%         mag = fftshift(abs(fft((rx_sig).^order,NFFT)));
%         [val, idx] = max(mag);
% 
%         offset = ff(idx) / order;
% 
%         figure; plot(ff,mag);
% %         figure; plot(real(rx_sig(1:114000)),'r')
% %         hold on; plot(imag(rx_sig(1:114000)),'b');
%         len = length(rx_sig);
%         rx_sig = rx_sig .* exp(-1i*2*pi*(offset)*(0:len-1)/sam_freq);
     
if isDQPSK
    %% Differential decoding
    Phase_Rotation = pi/4;
    % Normalization factor to convert from PI domain to linear domain
    Norm_Factor = modulationOrder/(2*pi);

    % Calculate the phase difference
    zPi = diff(unwrap([zeros(1, 1), angle(rx_sig)])) - Phase_Rotation;

    % Convert zPi to linear domain; and make hard decisions by rounding to the
    % nearest integer
    % Note: To be consistent with the blocks, we map 0.5 to 0
    demodSym = ceil(zPi*Norm_Factor - 0.5);

    % Remap to 0:M-1
    demodSym(demodSym < 0) = modulationOrder + demodSym(demodSym < 0);
    
    mappingTable =[0 0 1 1;...
                   0 1 0 1];
%         mappingTable =[0 1 0 1;...
%                        0 0 1 1];
    bitStream = mappingTable(:,demodSym+1 );
    demodSig = bitStream(:);
else
    %% carrier phase recovery, 1 sample per symbol
    %% incorrect method
    
%     numSamples = 100;
%     %         hPhaseSync = comm.PSKCarrierPhaseSynchronizer(modulationOrder,...
%     %                               'ObservationInterval',numSamples);
%     nSym = floor( length(rx_sig)/numSamples );
%     rx_sig = rx_sig(:);
%     phaseRcv = zeros(1,nSym*numSamples);
%     for i = 1:nSym
%         idx = (i-1)*numSamples+1: i*numSamples;
%         %       scatterplot(rx_sig(idx));
%         %               [rx_sig(idx), phEst(i)] = step(hPhaseSync, rx_sig(idx));   
%         phEst2 = 1/modulationOrder * angle(sum(rx_sig(idx).^modulationOrder)); % in degree
%               scatterplot(rx_sig(idx));
%         %       figure; plot(real(rx_sig(idx)),'-r.'); hold on; plot(imag(rx_sig(idx)),'-b.');
%         phaseRcv(idx) = rx_sig(idx).*exp(-1i*phEst2);
%          scatterplot(phaseRcv(idx));
%     end
    

    %% Demodulation
    rx_sig = rx_sig .* exp(-1i*initialPhase); % Dose AMR estimate initial phase?
    scatterplot(rx_sig);
     
    normFactor = modulationOrder/(2*pi); % normalization factor to convert from PI-domain to
                           % linear domain
    % convert input signal angle to linear domain; round the value to get ideal
    % constellation points 
    demodSym = round((angle(rx_sig) .* normFactor));

    % move all the negative integers by Modulation order
    demodSym(demodSym < 0) = modulationOrder + demodSym(demodSym < 0);  

    % convert symbols to bits
    %         demodSym = rx_sig;
    switch modulationOrder
    case 2
        demodSig = demodSym;
    case 4
        mappingTable =[0 0 1 1;...
                       0 1 0 1];
%         mappingTable =[0 1 0 1;...
%                        0 0 1 1];
        bitStream = mappingTable(:,demodSym+1 );
        demodSig = bitStream(:);
    case 8
        mappingTable =[0 0 0 0 1 1 1 1 ; ...
                       0 0 1 1 0 0 1 1 ; ...
                       0 1 0 1 0 1 0 1];
        bitStream = mappingTable(:,demodSym+1 );
        demodSig = bitStream(:);       
    end
end


        
catch err
    rethrow(err);  
end

function y = intanddump(x, Nsamp)
%INTANDDUMP Integrate and dump.
%   Y = INTANDDUMP(X, NSAMP) integrates the signal X for 1 symbol period, then
%   outputs one value into Y. NSAMP is the number of samples per symbol.
%   For two-dimensional signals, the function treats each column as 1
%   channel.
%
% --- Assure that X, if one dimensional, has the correct orientation --- %
wid = size(x,1);
if(wid ==1)
    x = x(:);
end

[xRow, xCol] = size(x);
x = mean(reshape(x, Nsamp, xRow*xCol/Nsamp), 1);
y = reshape(x, xRow/Nsamp, xCol);      

% --- restore the output signal to the original orientation --- %
if(wid == 1)
    y = y.';
end

function pushbutton_play_Callback(hObject, eventdata, handles)
sam_freq              = handles.sam_freq;
rx_sig                = handles.rx_sig;
decoder               = handles.decoder;
modType               = handles.modType;
soundFs               = handles.soundFs;
switch modType
    case {'AM','FM'}
        sound_out = rx_sig;
    otherwise
        
%% cvsd decoding
    switch decoder
        case 'CVSD'
            sw_enc_bit_dly_1 = 0;
            sw_enc_bit_dly_2 = 0;
            sw_step = 0;
            sw_rcnstrct_dly_1 = 0;
            sw_rcnstrct_dly_2 = 0;
            for i = 1: length(rx_sig)
                sw_enc_bit = rx_sig(i);
                [cvsd_decode_data,sw_enc_bit_dly_1,sw_enc_bit_dly_2,...
                    sw_step,sw_rcnstrct_dly_1,sw_rcnstrct_dly_2] = ...
                    cvsd_decode(sw_enc_bit,sw_enc_bit_dly_1,sw_enc_bit_dly_2,sw_step,sw_rcnstrct_dly_1,sw_rcnstrct_dly_2);
                sound_out(i) = cvsd_decode_data;
            end
    end
end


%% Play sound
try
    sound(sound_out,soundFs);
catch
    msgbox('Invalid sampling frequency! Control the resampling.');
end

%% update GUI states
set(handles.edit_viewStartIdx,'String','1');
set(handles.edit_viewEndIdx,'String',num2str(length(sound_out)));

%% save handle structure
handles.rx_sig          = sound_out;
handles.viewEndIdx      = length(sound_out);
handles.viewStartIdx    = 1;
pushbutton_viewWaveform_Callback(hObject, eventdata, handles);
     
function pushbutton_return_Callback(hObject, eventdata, handles)
handles.rx_sig = handles.defaultRx_sig;

%% update GUI states
updateLoggingMsg('Return proessed signal to original signal',handles);

%% save handle structure
guidata(hObject,handles);

function pushbutton_viewWaveform_Callback(hObject, eventdata, handles)

viewEndIdx      = handles.viewEndIdx;
viewStartIdx    = handles.viewStartIdx;
rx_sig          = handles.rx_sig;
sam_freq        = handles.sam_freq;

sam_period      = 1 / sam_freq;
displayRange    = viewStartIdx : viewEndIdx; %length(rx_sig); %handles.displayRange;

%% plot waveform
cla(handles.axes_demodMonitoring,'reset');
legend(handles.axes_demodMonitoring,'off');
hold(handles.axes_demodMonitoring, 'on');

tt = (0 : length(displayRange)-1)*sam_period ; 
plot(tt,real(rx_sig(displayRange)),'-b.','parent',handles.axes_demodMonitoring);
if ~isreal(rx_sig)
    plot(tt,imag(rx_sig(displayRange)),'-r.','parent',handles.axes_demodMonitoring);
    legend({'Inphase','Quadrature'});
end
xlim(handles.axes_demodMonitoring,[tt(1) tt(end)]);
xlabel(handles.axes_demodMonitoring,'Time[s]');
grid(handles.axes_demodMonitoring,'on');

function pushbutton_eyePattern_Callback(hObject, eventdata, handles)
rx_sig       = handles.rx_sig;
sps          = handles.samplesPerSymbol;
viewStartIdx = handles.viewStartIdx;
viewEndIdx   = handles.viewEndIdx;

displayRange    = viewStartIdx : viewEndIdx; %length(rx_sig); %handles.displayRange;
eyediagram(real(rx_sig(displayRange)),sps);

function pushbutton_scatterPlot_Callback(hObject, eventdata, handles)
rx_sig       = handles.rx_sig;
viewStartIdx = handles.viewStartIdx;
viewEndIdx   = handles.viewEndIdx;
displayRange    = viewStartIdx : viewEndIdx;

scatterplot(rx_sig(displayRange));

function pushbutton_viewSpectrum_Callback(hObject, eventdata, handles)

rx_sig         = handles.rx_sig;      
sam_freq       = handles.sam_freq;
viewEndIdx     = handles.viewEndIdx;
viewStartIdx   = handles.viewStartIdx;

rx_sig = rx_sig(viewStartIdx : viewEndIdx );
sigLen = length(rx_sig);
rx_sig = rx_sig(:);

FFTlength = 2.^nextpow2(sigLen);
%% zero-padding
if sigLen < FFTlength
    rx_sig = [rx_sig; zeros(FFTlength-(viewEndIdx-viewStartIdx+1) ,1) ];
    %set(handles.text_monitor_msg,'String','Signal is zero-padded');
end

%% Calculated power spectral density
% [Pxx,~] = pwelch(rx_sig,window,noverlap,FFTlength, ...
%     sam_freq,'twosided');
window           = hanning( FFTlength );              % Time domain data tapering window
noverlap         = FFTlength/2;                     % The number of samples for overlap
freqResoultion   = sam_freq / FFTlength;

numberOfSegments = (length(rx_sig)-noverlap)./(FFTlength-noverlap);
numberOfSegments = fix(numberOfSegments);

LminusOverlap = FFTlength-noverlap;
xStart = 1:LminusOverlap:numberOfSegments*LminusOverlap;
xEnd   = xStart+FFTlength-1;

% accumulatedSxx = zeros(numberOfSegments, FFTlength);
Sxx = 0; Sxxk=0;

for i = 1:numberOfSegments            
    xw = rx_sig(xStart(i):xEnd(i)).*window; % Window the data
    U = window'*window;                     % compensates for the power of the window.
    Xx = fft(xw,FFTlength);
    Sxxk = Xx.*conj(Xx)/U;                  % Auto spectrum.
    Sxx  = Sxx + Sxxk;
end     
Pxx = Sxx;

% To get 0 dB reference point 
% Pxx = Pxx./max(Pxx);
Pxx = fftshift(Pxx);          Pxx_dB =( 10*log10(Pxx) ); 

% set frequency axis
freqAxis = ( -FFTlength/2 : (FFTlength -1)/2 ) / FFTlength * sam_freq;

%% plot spectrum
cla(handles.axes_demodMonitoring,'reset');
legend(handles.axes_demodMonitoring,'off');
hold(handles.axes_demodMonitoring, 'on');

plot(freqAxis, Pxx_dB,'parent',handles.axes_demodMonitoring);
xlabel('Frequency (Hz)','parent',handles.axes_demodMonitoring);
ylabel('Magnitude (dB/Hz)','parent',handles.axes_demodMonitoring);
xlim(handles.axes_demodMonitoring,[freqAxis(1), freqAxis(end)]);
grid(handles.axes_demodMonitoring,'on'); 

%% Update handles structure
handles.freqAxis = freqAxis;
handles.Pxx = Pxx;

guidata(hObject,handles);

function edit_viewStartIdx_Callback(hObject, eventdata, handles)
viewStartIdx = str2double( get(hObject,'String') );
if viewStartIdx <= 1
    viewStartIdx = 1;
end
set(hObject,'String', num2str(viewStartIdx));
handles.viewStartIdx = viewStartIdx;
guidata(hObject,handles);

function edit_viewStartIdx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_viewEndIdx_Callback(hObject, eventdata, handles)
rx_sig       = handles.rx_sig;
viewEndIdx   = str2double( get(hObject,'String') );
sigMaxLength = length(rx_sig);

if viewEndIdx >= sigMaxLength
    viewEndIdx = sigMaxLength;
end

set(hObject,'String', num2str(viewEndIdx));

handles.viewEndIdx = viewEndIdx;
guidata(hObject,handles);

function edit_viewEndIdx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_numberOfInputSamples_Callback(hObject, eventdata, handles)

function edit_numberOfInputSamples_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = popupmenu_resampleRate_Callback(hObject, eventdata, handles)

viewStartIdx = handles.viewStartIdx;
viewEndIdx   = handles.viewEndIdx;
rx_sig       = handles.rx_sig;
sam_freq     = handles.sam_freq;

sigLen = length(rx_sig);
%% get data from popupmenu
contents = cellstr(get(hObject,'String'));
resampleRate = contents{get(hObject,'Value')};

slashIdx = strfind(resampleRate,'/');

if isempty(slashIdx)
    p =str2double(resampleRate);
    q = 1;    
    if ~(p == 1 && q == 1)
        rx_sig = interp(rx_sig,p);
    end
else
    p = str2double(resampleRate(1:slashIdx-1));
    q = str2double(resampleRate(slashIdx+1:end));
    if ~(p == 1 && q == 1)
        rx_sig = decimate(rx_sig,q);
    end
    
end
myMsg = ['Change sampling rate by ',num2str(p),'/',num2str(q),' times'];
updateLoggingMsg(myMsg,handles);
% filename(startIdx:end) = [];
%% resample the data by p/q times
% rx_sig = resample(rx_sig,p,q);
viewEndIdx = ceil(sigLen*p/q);

sam_freq =  sam_freq *p/q;
%% update handles structure
handles.resampleRate        = resampleRate;
handles.rx_sig              = rx_sig;
handles.viewStartIdx        = viewStartIdx;
handles.viewEndIdx          = viewEndIdx;
handles.sam_freq            = sam_freq;
guidata(hObject,handles);
  
%% update GUI states
set(handles.edit_viewStartIdx,'String',viewStartIdx);
set(handles.edit_viewEndIdx,'String',viewEndIdx);
set(handles.edit_numberOfInputSamples,'String',length(rx_sig));
set(handles.edit_samFreq,'String',num2str(sam_freq));
cla(handles.axes_demodMonitoring,'reset');

function popupmenu_resampleRate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_modType_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String')); %returns popupmenu_modType contents as cell array
modType = contents{get(hObject,'Value')}; %returns selected item from popupmenu_modType

%% update GUI states
myMsg =['Selected demodulation type is ', modType];
updateLoggingMsg(myMsg,handles);

%% save handles structure
handles.modType = modType;
guidata(hObject,handles);

function popupmenu_modType_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% update monitoring message box
function updateLoggingMsg(myMsg,handles)
oldmsgs = cellstr(get(handles.text_demod_monitor_msg,'String'));
if length(oldmsgs) >= 8
    set(handles.text_demod_monitor_msg,'String',...
        [   oldmsgs{end-6}; ...
            oldmsgs{end-5}; ...
            oldmsgs{end-4}; ...
            oldmsgs{end-3}; ...
            oldmsgs{end-2}; ...
            oldmsgs{end-1}; ...
            oldmsgs{end}; ...
            {myMsg}...
         ] );
%                  oldmsgs{end}; ...
else
    set(handles.text_demod_monitor_msg,'String',[oldmsgs;{myMsg}] );
end
pause(0.0000001);       

function pushbutton_carrierFreqOffset_Callback(hObject, eventdata, handles)
rx_sig     = handles.rx_sig;
sam_freq   = handles.sam_freq;
sam_period = 1 / sam_freq;
coarseFreqOffset = handles.coarseFreqOffset;

ttIdx = ( 0 : length(rx_sig)-1) * sam_period;
ttIdx = ttIdx(:);
rx_sig = rx_sig(:);
rx_sig = rx_sig.*exp(-1i*(2*pi*coarseFreqOffset*ttIdx ) );

handles.rx_sig = rx_sig;
guidata(hObject,handles);

function pushbutton_reset_Callback(hObject, eventdata, handles)
close(gcbf);    Demodulation;

function edit_sps_Callback(hObject, eventdata, handles)
sps   = str2double( get(hObject,'String') );
handles.samplesPerSymbol = sps;
guidata(hObject,handles);

function edit_sps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_soundFs_Callback(hObject, eventdata, handles)
soundFs   = str2double( get(hObject,'String') );
handles.soundFs = soundFs;
guidata(hObject,handles);
updateLoggingMsg(['Set sampling frequency for sound play : ',num2str(handles.soundFs),' Hz'],handles);

function edit_soundFs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [cvsd_decode_data,sw_enc_bit_dly_1,sw_enc_bit_dly_2,sw_step,sw_rcnstrct_dly_1,sw_rcnstrct_dly_2] = cvsd_decode(sw_enc_bit,sw_enc_bit_dly_1,sw_enc_bit_dly_2,sw_step,sw_rcnstrct_dly_1,sw_rcnstrct_dly_2)

% 07.11.5
% lidongshi
% cvsd decode
%
DELTA_MIN = 0.0002;%4;		% 0.0002 Scaled to 2 ^ 14
DELTA_MAX = 0.0078;%128;     % 0.0078 Scaled to 2 ^ 14
SYLLABIC_CONST = 0.9845;%32260;   % 0.9845 Scaled to 2 ^ 15
INTEG_B1 = 1.2708;%20821;   % 1.2708 Scaled to 2 ^ 14
INTEG_B2 = 0.3202;%10492;   % 0.3202 Scaled to 2 ^ 15
INTEG_G2D = 1.5092;%24726;   % 1.5092 Scaled to 2 ^ 14

% sw_enc_bit_dly_1 = 0;
% sw_enc_bit_dly_2 = 0;
% sw_step = 0;
% sw_rcnstrct_dly_1 = 0;
% sw_rcnstrct_dly_2 = 0;

% 	int sw_acc;
% 	int sw_rcnstrct;
% 	
% 	long int tmp1;
% 	long int tmp2;
% 	long int tmp3;
% 	int tmp4;
tmp_lds_1 = 0;
tmp_lds_2 = 0;
tmp_lds_3 = 0;
%  Judging continuous three bit identity and syllablic filter
sw_enc_bit = 2 * sw_enc_bit - 1;   % 0 --> -1; 1 --> 1

sw_acc = sw_enc_bit + sw_enc_bit_dly_1 + sw_enc_bit_dly_2;

tmp_lds_1 = sw_step;
tmp1 = SYLLABIC_CONST * tmp_lds_1;
% 	//tmp1 = (SYLLABIC_CONST * sw_step);
% 	tmp1 >>= 15;
	
if abs(sw_acc) == 3
	sw_step = tmp1 + DELTA_MAX;
else
	sw_step = tmp1 + DELTA_MIN;
end
% 	Primary reconstruction integration
tmp_lds_2 = sw_rcnstrct_dly_1;
tmp1 = INTEG_B1 * tmp_lds_2;
% 	//tmp1 = INTEG_B1 * sw_rcnstrct_dly_1;
% 	tmp1 >>= 14;
tmp_lds_2 = sw_rcnstrct_dly_2;
tmp2 = INTEG_B2 * tmp_lds_2;
% 	//tmp2 = INTEG_B2 * sw_rcnstrct_dly_2;
% 	tmp2 >>= 15;
tmp_lds_1 = sw_step;
tmp3 = INTEG_G2D * tmp_lds_1;
% 	//tmp3 = INTEG_G2D * sw_step;
% 	tmp4 = tmp3 >> (14 + 1);           % '14' add '1' for normalization
tmp3 = tmp3 * sw_enc_bit;

tmp1 = tmp1 - tmp2 + tmp3;

% 	Saturation process
if (tmp1 >= 32767)        % up overflow
	sw_rcnstrct = 32767;
elseif (tmp1 <= -32768)  % down overflow
	sw_rcnstrct = -32768;
else                      % no overflow
	sw_rcnstrct = tmp1;
end
	
% 	Shift register shift
sw_enc_bit_dly_2 = sw_enc_bit_dly_1;
sw_enc_bit_dly_1 = sw_enc_bit;

sw_rcnstrct_dly_2 = sw_rcnstrct_dly_1;
sw_rcnstrct_dly_1 = sw_rcnstrct;
	
% return (sw_rcnstrct);
cvsd_decode_data = sw_rcnstrct;
