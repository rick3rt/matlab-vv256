% RECONSTRUCTION OF HIGH FRAME RATE IMAGING
% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% Description: 
%   Sequence programming file for L12-5_38mm Linear array, using 2-1 synthetic   
%   aperture plane wave transmits and receive acquisitions on 128 channels 
%   system. 128 transmit channels and 96 receive channels are active and 
%   positioned as follows (each char represents 4 elements) for each of the 
%   2 synthetic apertures.
%
%   Element Nos.                                1              1
%                               6       9       2              9
%               1               5       7       9              2
%   Aperture 1: |               |       |       |              |
%               tttttttttttttttttttttttttttttttt----------------
%               rrrrrrrrrrrrrrrrrrrrrrrr------------------------
%   Aperture 2: |               |       |       |              |
%               ----------------tttttttttttttttttttttttttttttttt
%               ------------------------rrrrrrrrrrrrrrrrrrrrrrrr
%               |               |       |       |              | 
%
%   The receive data from each of these apertures are stored under  
%   different acqNums in the Receive buffer. The reconstruction sums the 
%   IQ data from the 2 aquisitions and computes intensity values to produce 
%   the full frame. Processing is asynchronous with respect to acquisition.
%
% Last update 
%    12/13/15   update to SW 3.0
%    08/01/19   modifications by Rick W.

clear; clc; close all;
% suppress matlab warnings:
%#ok<*SAGROW>   % variable size changes in loop
%#ok<*UNRCH>    % cannot reach code after return


%% load user parameters
P = load_parameters('IQR');
disp(P)

%% GENERAL SETTINGS AND PARAMETERS
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.fakeScanhead = 1; % optional (if probe not conneted)
Resource.Parameters.simulateMode = 2; % force always mode 2
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.


%% TRANSDUCER
% Specify Trans structure array.
Trans.name = 'L12-5 38mm';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L12-5_50mm transducer is 'known' transducer so we can use computeTrans.


%% PDATA IMAGE RECONSTRUCTION (REAL TIME)
% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.


%% MEDIA (SIMULATION)
% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';


%% RESOURCE BUFFER
% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*2; % this size allows for 3 acqs, maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numOfFrames;        % 40 frames used for RF cineloop.

Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = P.numOfFrames; % IQDATA IS HERE
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;

Resource.DisplayWindow(1).Title = 'L12-5_38mmFlashs';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = P.numOfFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);


%% TRANSMIT WAVEFORM
% Specify Transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,1,1];   % A, B, C, D  for 8.1 MHz transmit frequency


%% TRANSMIT EVENTS
% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, 3);
TX(2).aperture = 65;  % Use the tx aperture that starts at element 65.


%% TIME GAIN CONTROL
% Specify TGC Waveform structure.
TGC.CntrlPts = [234,368,514,609,750,856,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);


%% RECEIVE STRUCTURE
% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,128), ...
                        'aperture', 1, ...
                        'startDepth', P.startDepth, ...
                        'endDepth',maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,2*Resource.RcvBuffer(1).numFrames);
%

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames  % 2 acquisitions per frame
    Receive(2*i).callMediaFunc = 1;
    % -- 1st synthetic aperture acquisition for aperture 1.
    Receive(2*i-1).Apod(1:96) = 1.0;
    Receive(2*i-1).aperture = 1;
    Receive(2*i-1).framenum = i;
    Receive(2*i-1).acqNum = 1;
    Receive(2*i-1).callMediaFunc = 1;
    % -- 2nd synthetic aperture acquisition for aperture 65.
    Receive(2*i).Apod(33:128) = 1.0;
    Receive(2*i).aperture = 65;
    Receive(2*i).framenum = i;
    Receive(2*i).acqNum = 2;
end


%% RECONSTUCTION SETTINGS
% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:2),1,Resource.RcvBuffer(1).numFrames);

% make reconstrucion event per frame
for i = 1:Resource.RcvBuffer(1).numFrames
    Recon(i).IntBufDest= [1,i];         % each frame in its own interbuffer
    Recon(i).RINums = (i-1)*2 + (1:2);
    Recon(i).rcvBufFrame = i;
end

% Define ReconInfo structures (depent on frame number).
for i = 1:Resource.RcvBuffer(1).numFrames
    ReconInfo((i-1)*2+(1:2)) = repmat(struct('mode', 'replaceIQ', ...  % replace IQ data.
           'txnum', 1, ...
           'rcvnum', (i-1)*2 +1, ...
           'regionnum', 1), 1, 2);

    % - Set specific ReconInfo attributes.
    ReconInfo((i-1)*2 +2).mode = 'accumIQ_replaceIntensity'; % accumulate and detect IQ data in output buffer.
    ReconInfo((i-1)*2 +2).txnum = 2;
    ReconInfo((i-1)*2 +2).rcvnum = (i-1)*2 +2;
end


%% PROCESSING OF THE DATA (SHOW IMAGE)
% Specify Process structure array.
pers = P.persistLevel; %20
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',6.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
%

%% EVENTS AND SEQUENCE CONTROL
% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'pause';
SeqControl(2).condition = 'extTrigger'; % input BNC #1 falling edge
SeqControl(2).argument = 1; % 500 msec timeout delay
    % (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
SeqControl(3).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(3).argument = 200;%200;  % 200 usec
SeqControl(4).command = 'timeToNextAcq';  % time between frames
SeqControl(4).argument = 9000 - 200;  % 20000 usec = 20 msec
SeqControl(5).command = 'returnToMatlab';

nsc = 6; % nsc is count of SeqControl objects
n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = '1st aperture.';
    Event(n).tx = 1;        % use TX structure of 1st aperture.
    Event(n).rcv = 2*i-1;    % use Rcv structure of 1st aperture.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 3; % time between syn. aper. acqs.
    n = n+1;
    
    Event(n).info = '2nd aperture.';
    Event(n).tx = 2;       % use Tx structure of 2nd aperture.
    Event(n).rcv = 2*i;    % use Rcv structure of 2nd aperture.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = nsc; % time between syn. aper. acqs.
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = i;      % reconstruction
    Event(n).process = 1;    % processing
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame 
        Event(n).seqControl = 5;
    end
    n = n+1;
end

% Event(n).info = 'Jump back to first event';
% Event(n).tx = 0;        % no TX
% Event(n).rcv = 0;       % no Rcv
% Event(n).recon = 0;     % no Recon
% Event(n).process = 0; 
% Event(n).seqControl = 1; % jump command


%% UI CONTROL
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
frameRateFactor = P.frameRateFactor; %2


%% SAVE AND RUN VSX
% Save all the structures to a .mat file.
save(['MatFiles/',filename]);

% LOAD RCV DATA TO RECONSTRUCT
disp('LOADING RCVDATA...')
load('testRcvData.mat')
disp('DONE LOADING RCVDATA!')
beep

% run VSX
VSX

%% CODE TO RUN AFTER VSX CLOSES
beep
answer = input('Save IQData? [1: yes, 0: no]');
if answer == 1
    save('testIQData.mat', 'IQData');
    disp('DONE SAVING!')
    beep
end

%% Place this return to prevent executing code below
return


%% DEFINITION OF FUNCTIONS (UI AND EF)
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
evalin('base','if VDAS==1, Result = loadTgcWaveform(1); end');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');

return
%RangeChangeCallback


% that's it