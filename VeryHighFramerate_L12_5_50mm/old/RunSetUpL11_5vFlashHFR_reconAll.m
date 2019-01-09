% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name RunSetUpL11_5vFlashHFR_reconAll.m: "High Frame Rate" post-acquisition processing script (script 2 of 2)
%
% Description: 
% To support higher acquisition frame rates with reduced DMA overhead, the first script (SetupL11_5vFlashHFR_acquire) acquires
%   numAcqs acquisitions into each RcvBuffer "super" frame, performing a transferToHost after each numAcqs set. After quitting VSX, 
%   the RF data is saved in a matfile called 'RFdataHFR'. The real-time image reconstruction only processes the first subframe  
%   in the super frame to minimize processing time and yet provide real-time feedback during data collection. 
%   The data is saved in a matfile to permit the "clear all" at the beginning of every script. Note that in this HFR script, 
%   each acquisition will lead to an image frame, and each DMA data transfer to the host contains one receive "super" frame.
% 
% This second script has been written to reconstruct and displays all acquisitions. To be compatible with the acquired data, 
%   it must use the same RcvBuffer size and is intended to be run in simulateMode=2 using the previously collected data. It works
%   in real-time as well, but the processing cannot keep up with acquisition, and some superframes are skipped but the script 
%   displays all sub-frames within a super frame. 
% 
% To illustrate the HFR effect in simulation, the 'movePoints' function is invoked (in SetupL11_5vFlashHFR_acquire) 
%   for every acquisition, resulting in a very discontinous apparent motion of the scatterers in "real-time" when processing 
%   only one acquisition per superframe. However, the motion is smooth on playback when processing each T/R acquisition.
%
% last modification 12/13/2015

clear all

simulateMode = 2; % use mode=2 for post processing data collected by 'SetupL11_5vFlashHFR' and stored in 'RFdataHFR.mat'.
                  % use mode=0 for acquisiton, and processing of every acquisition until superframes are skipped

if simulateMode ==2,
    load ('RFdataHFR')  % this file includes RcvData and the numAcqs and numFrames parameters for properly setting Resources
else
    P.startDepth = 5;   % Acquisition depth in wavelengths
    P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.
    P.numAcqs = 100;  % default no. of acquisitions in a frame (this is a "superframe") for testing
    P.numFrames = 10;  % default no. of frames    
end

% Define system parameters.
filename = mfilename;
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numAcqs;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;       % one frame is defined by the amount of data in 1 'transfer to host'
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = P.numFrames*P.numAcqs; % make a cine loop that can hold images from the entire data buffer
Resource.DisplayWindow(1).Title = filename;
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = P.numFrames*P.numAcqs;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];   

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,298,395,489,618,727,921,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));

Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1),1,P.numFrames*P.numAcqs);
                    
% - Set event specific Receive attributes.
for frame = 1:P.numFrames    % -- Acquisitions for 'super' frame.
    Receive(P.numAcqs*(frame-1) + 1).callMediaFunc = 1;   % move media points only every super frame
    for acqNum = 1:P.numAcqs
        rcvnum = P.numAcqs*(frame-1) + acqNum;
        Receive(rcvnum).Apod(:)  = 1;           
        Receive(rcvnum).framenum = frame;       % each frame is transfered by one DMA ("transfer to host").
        Receive(rcvnum).acqNum   = acqNum;      % each acquisition is stored separately, indexed by framenum and acqNum
    end
end

% Specify Recon structure arrays: only need to define Recon for one Receive frame (superframe)
% - We need one Recon structure for each image to be produced, and here define all Recons needed for one superframe.
% In this Flash acquisition script, each acquisition is reconstructed and displayed as an image, and therefore
% each Recon will need just one ReconInfo structure, since every acquisition will be used to create a unique image.
% So, define ReconInfos in the same loop as Recons
for acqNum = 1:P.numAcqs
    Recon(acqNum) = struct('senscutoff', 0.6, ...
                    'pdatanum', 1, ...
                    'rcvBufFrame', -1, ...      % use most recently transferred frame (in simulateMode=2, the frames are processed seqeuntially)
                    'IntBufDest', [1,1], ...
                    'ImgBufDest', [1,-1], ...   % write to a different ImageBuffer each recon using autoincrement
                    'RINums', acqNum);          % one unique ReconInfo per Recon 
                                                % NOTE: can never use the same ReconInfo in different Recons!
    
    ReconInfo(acqNum) = struct('mode', 'replaceIntensity', ...    
                    'txnum', 1, ...
                    'rcvnum', acqNum, ...       % the index of the receive acquisition
                    'regionnum', 1);
end
               

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
%  for acqNum = 1:P.numAcqs                   % Use this loop if independent Process structures are desired for explicit framenum
%      Process(acqNum) = Process(1);        % designation. The specific process structures must be referenced in the event sequence. 
%      Process(acqNum).Parameters{4} = acqNum;
%  end
 
% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 200;  % 200 usecs
SeqControl(3).command = 'triggerOut';
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'noop';
SeqControl(5).argument = 50000; % 10 ms
nsc = 6; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for frame = 1:P.numFrames
    for acqNum = 1:P.numAcqs  % Acquire all acquisitions defined in one super frame
        Event(n).info = 'Acquire RF';
        Event(n).tx = 1;         
        Event(n).rcv = P.numAcqs*(frame-1) + acqNum;    
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = [2,3]; 
        n = n+1;
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [2,3,nsc];
        SeqControl(nsc).command = 'transferToHost'; % defines the receive frame 
        nsc = nsc + 1;
    
    % Do reconstruction and processing for each acq in super frame
    for acqNum = 1:P.numAcqs
        Event(n).info = 'Reconstruct';
        Event(n).tx = 0;         
        Event(n).rcv = 0;        
        Event(n).recon = acqNum; 
        Event(n).process = 1;    
        Event(n).seqControl = 0;
        n = n+1;
    end
    Event(n-1).seqControl = 4;
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        
Event(n).rcv = 0;       
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = 1; 


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

% Save all the structures to a .mat file.
save(['MatFiles/',filename]); 
disp([ mfilename ': NOTE -- Running VSX automatically!'])
VSX    

return


% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);
evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
