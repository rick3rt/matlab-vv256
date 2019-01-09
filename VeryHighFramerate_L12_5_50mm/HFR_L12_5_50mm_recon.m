%% HIGH FRAME RATE IQDATA RECONSTRUCTION
% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name: RunSetUpL11_5vFlashHFR_acquire.m - "High Frame Rate"
%   acquisition for ultrafast imaging (script 1 of 2).
%
% Description: NOTE: This is a unique example that runs VSX automatically,
%   and after closing the Control GUI, runs the post-processing script
%   (script 2 of 2) automatically. Two scripts are part of this example!!!
%
% To support higher acquisition frame rates with reduced DMA overhead, this script 
%   acquires a large set of T/R acquisitions into each RcvBuffer 'super'
%   frame, performing a transferToHost only after each group of "numAcqs"
%   acquisitions. Even though each acquisition is intended to be
%   reconstructed and processed into a display frame, live image
%   reconstruction only processes the first acquisition in the super frame
%   in order not to slow down real-time acquisition, and yet provide visual
%   feedback to help the user collect desired data. To post-process all
%   acquisitions and superframes present in the RcvBuffer, the data is
%   saved in a matfile when VSX quits (see end of the script).
%   Unfortunately, the current system software does not permit running this
%   same script to process all of the acquisitions, even in playback
%   (simulationMode=2), because a sequence can only reconstruct a maximum
%   of 16 acquisitions in a given receive frame, and therefore cannot
%   include all of the data transferred in one DMA (here called a super
%   frame).
% 
% A second script for reconstructing and displaying all aquisitions has
%   been developed to be run in simulateMode=2 after this one,
%   using the previously collected data. That script,
%   SetUpL11_5vFlashHFR_reconAll.m, reconstructs and displays all
%   acquisitions by defining ReconInfo structures for each acquisition, and
%   defining a Recon event for each display frame. Acquiring data using
%   that script results in skipped superframes when processing cannot keep
%   up with acquisition.
%
% To illustrate the HFR acquisition and display in simulation, the
%   'movePoints' function is invoked for every acquisition, resulting in a
%   very discontinous apparent motion of the scatterers in "real-time".
%   However, the motion is smooth on playback when post-processing each T/R
%   acquisition with the "reconAll" script.
%
% For convenience, this script is currently set to launch VSX
%   automatically. Upon quitting VSX (by closing the control GUI), the
%   "reconAll" script is automatically invoked to process and display all
%   of the frames just collected. Note that in simulateMode=1, this script
%   runs once, filling all numFrames frames, then freezing. Simply close
%   the GUI, and the reconAll script launches automatically to demonstrate
%   the post-processing of all acquisitions.
%
% Last modification 12/13/2015
% 
% RW: 09/01/2019 rewrite for L12-5 50mm

clear; clc; close all;
% suppress matlab warnings:
%#ok<*SAGROW>   % variable size changes in loop
%#ok<*UNRCH>    % cannot reach code after return


%% load user parameters
P = HFR_load_parameters('recon');


%% GENERAL SETTINGS AND PARAMETERS
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = P.simMode;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.


%% TRANSDUCER
% Specify Trans structure array.
Trans.name = 'L12-5 50mm';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % L12-5 50mm transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;      % set maximum high voltage limit for pulser supply.


%% PDATA IMAGE RECONSTRUCTION (REAL TIME)
% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.


%% MEDIA (SIMULATION)
% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';


%% RESOURCE BUFFER
% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numAcqs; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numSuperFrames;        % Number of 'Super Frames'

Resource.InterBuffer(1).numFrames = P.numSuperFrames*P.numFramesPerSuperFrame; % will contain IQData
% Resource.InterBuffer(1).datatype = 'complex';

Resource.ImageBuffer(1).numFrames = 10;
% Resource.ImageBuffer(1).datatype = 'double';

Resource.DisplayWindow(1).Title = 'L12-5_50mmFlash HFR';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
% Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = P.numSuperFrames*P.numAcqs; %P.numSuperFrames; (?)
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);


%% TRANSMIT WAVEFORM
% Specify Transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,1,1];   % A, B, C, D  for 8.1 MHz transmit frequency
% TW(1).Parameters = [Trans.frequency,.67,2,1];   % A, B, C, D  for 8.1 MHz transmit frequency

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
TX(3).aperture = 129; % Use the tx aperture that starts at element 129.


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
                        'callMediaFunc', 1),1,P.numSuperFrames*P.numAcqs);
% movepoints EVERY acquisition to illustrate superframe concept
% real-time images will look "jerky" but using the reconstructAll script,
% playback process all acquisitions and shows smooth displacement

% Apodization factor per acquisition (3 apertures)
ApodFactors = zeros(3,128);
ApodFactors(1,1:85) = 1;
ApodFactors(2,22:107) = 1;
ApodFactors(3,44:128) = 1;

% Aperture HVMux selection (3 apertures)
ApertureList = [1,65,129];

% - Set event specific Receive attributes.
for i = 1:P.numSuperFrames                  % number of Super Frames = P.numSuperFrames
    for j = 1:P.numFramesPerSuperFrame      % number of frames per superframe
        for k = 1:P.numApertures            % number of acquisitions (i.e. apertures) per frame
            % -- Acquisitions for 'super' frame.
            rcvNum = P.numFramesPerSuperFrame*P.numApertures*(i-1) + P.numApertures*(j-1) + k;
            % set Reveice objects
            Receive(rcvNum).Apod = ApodFactors(k,:); 
            Receive(rcvNum).aperture = ApertureList(k); 
            Receive(rcvNum).framenum = i; % superframe num
            Receive(rcvNum).acqNum = P.numApertures*(j-1)+k; % acq num in superframe
            % dont call media func in sim mode between apertures.
            if k ~= 3
                Receive(rcvNum).callMediaFunc = 0;
            end
        end
    end
end

% DEBUGGING HELP
% fprintf('% 5d\t% 5d\t% 5d\t% 5d\t% 5d\n', rcvNum, Receive(rcvNum).acqNum, i, j, k)


%% RECONSTUCTION SETTINGS
% Specify Recon structure arrays
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:3),1,P.numSuperFrames*P.numFramesPerSuperFrame);

% make reconstrucion event per frame
for i = 1:P.numSuperFrames                  % superFrames
    for j = 1:P.numFramesPerSuperFrame      % frames within superFrame
        reconNum = (i-1)*P.numFramesPerSuperFrame + j;
        Recon(reconNum).IntBufDest= [1, reconNum]; % each frame in its own interbuffer
        Recon(reconNum).RINums = (i-1)*P.numFramesPerSuperFrame + (j-1)*P.numApertures + (1:3);
        Recon(reconNum).rcvBufFrame = i;
    end
end

% initialize ReconInfo structures (depent on frame number).
ReconInfo = repmat(struct('mode', 'replaceIQ', ...  % replace IQ data.
               'txnum', 1, ...
               'rcvnum', 1, ...
               'scaleFactor', 2.0, ...
               'regionnum', 1), 1, P.numSuperFrames*P.numAcqs); % reconInfo per receive object
%

% reDefine ReconInfo structures
for i = 1:P.numSuperFrames                  % superFrames
    for j = 1:P.numFramesPerSuperFrame      % frames within superFrame
        % set ReconInfo per aperture
        % Aperture 1
        riNum = P.numFramesPerSuperFrame*P.numApertures*(i-1) + P.numApertures*(j-1) + 1;
        % fprintf('% 5d\n',riNum)
        ReconInfo(riNum).mode = 'replaceIQ';
        ReconInfo(riNum).txnum = 1;
        ReconInfo(riNum).rcvnum = riNum;
        
        % Aperture 2
        riNum = riNum + 1;
        ReconInfo(riNum).mode = 'accumIQ';
        ReconInfo(riNum).txnum = 2;
        ReconInfo(riNum).rcvnum = riNum;
        
        % Aperture 3
        riNum = riNum + 1;
        ReconInfo(riNum).mode = 'accumIQ_replaceIntensity';
        ReconInfo(riNum).txnum = 3;
        ReconInfo(riNum).rcvnum = riNum;
    end
end


%% PROCESSING OF THE DATA (SHOW IMAGE)
% Specify Process structure array.
pers = P.persistLevel; %20
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',4.0,...     % pgain is image processing gain
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
SeqControl(2).command = 'timeToNextAcq';  % time between Aperture acquisitions
SeqControl(2).argument = P.time.acq;  % 200 usecs
SeqControl(3).command = 'timeToNextAcq';  % time between individual Frame acquisitions
SeqControl(3).argument = P.time.frame;  % 200 usecs
SeqControl(4).command = 'triggerOut';
SeqControl(5).command = 'returnToMatlab';
SeqControl(6).command = 'timeToNextAcq';  % time to acquire new superframes
SeqControl(6).argument = P.time.superFrame;  % 10 msecs

nsc = 7;    % nsc is count of SeqControl objects
n = 1;      % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:P.numSuperFrames                  % number of Super Frames = P.numSuperFrames
    for j = 1:P.numFramesPerSuperFrame      % number of frames per superframe
        for k = 1:P.numApertures            % number of acquisitions (i.e. apertures) per frame
            % 3 apertures
            Event(n).info = ['Acquire RF - aperture ' num2str(k) ' - frame ' num2str(j) ' - superFrame ' num2str(i)] ;
            Event(n).tx = 1;         
            Event(n).rcv = P.numFramesPerSuperFrame*P.numApertures*(i-1) + P.numApertures*(j-1) + k; % same as rcvNum
            Event(n).recon = 0;      
            Event(n).process = 0;    
            Event(n).seqControl = 2; % timeToNextAcq - aperture. Only trigger when complete frame is recorded, not per aperture
            n = n+1;
        end
        % when the apertures are received, next frame within superFrame can have different timeToNectAcq
        Event(n-1).seqControl = [3,4]; % timeToNextAcq - frame. Also triggerOut after all apertures received
    end
    
    % When complete superFrame is recorded, transferToHost
    Event(n-1).seqControl = [6,4,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
    nsc = nsc + 1;
    
    % Do reconstruction and processing for each acq in super frame
    for j = 1:P.numFramesPerSuperFrame
        Event(n).info = ['Reconstruct - frame ' num2str(j) ' - superFrame ' num2str(i)]; 
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = P.numFramesPerSuperFrame*P.numApertures*(i-1) + j;
        Event(n).process = 1;
        Event(n).seqControl = 0; % seqControl = 0 means no operation (noop)
        n = n+1;
    end

end

% DEBUGGING
% fprintf('\n\n');
% for i = 1:numel(Event)
%     fprintf([Event(i).info '\n'])
% end


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
frameRateFactor = P.numAcqs;


%% SAVE AND RUN VSX
% Save all the structures to a .mat file.
save(['MatFiles/' filename]);

disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX    
commandwindow  % just makes the Command window active to show printout


%% CODE TO RUN AFTER VSX CLOSES


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
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback


% that's it