% File name SetUpL11_4vAcquireRF.m
% - Synchronous acquisition into a single RcvBuffer frame.
clear


%% GENERAL SETTINGS AND PARAMETERS
% Specify system parameters
Resource.Parameters.numTransmit = 128; % no. of transmit channels
Resource.Parameters.numRcvChannels = 128; % change to 64 for Vantage 64 system
Resource.Parameters.connector = 1; % trans. connector to use (V 256).
Resource.Parameters.speedOfSound = 1540; % speed of sound in m/sec
Resource.Parameters.fakeScanhead = 1; % optional (if no L7-4)

% simulate? 0 real thing, 1 sim mode
Resource.Parameters.simulateMode = 1; % runs script with hardware


%% MEDIA (SIMULATION)
% Specify media points
Media.MP(1,:) = [0,0,100,1.0]; % [x, y, z, reflectivity]


%% TRANSDUCER
% Specify Trans structure array.
Trans.name = 'L11-4v';
Trans.frequency = 6.25; % not needed if using default center frequency
Trans.units = 'wavelengths';
Trans = computeTrans(Trans); % L11-4v transducer is 'known' transducer.


%% RESOURCE BUFFER
% Specify Resource buffers.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048;
Resource.RcvBuffer(1).colsPerFrame = 128;
Resource.RcvBuffer(1).numFrames = 1; % minimum size is 1 frame.


%% TRANSMIT WAVEFORM
% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [5.208,0.67,2,1]; % A, B, C, D

% Specify TX structure array.
TX(1).waveform = 1; % use 1st TW structure.
TX(1).focus = 0;
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));


%% TIME GAIN CONTROL
% Specify TGC Waveform structure.
TGC(1).CntrlPts = [500,590,650,710,770,830,890,950];
TGC(1).rangeMax = 200;
TGC(1).Waveform = computeTGCWaveform(TGC);


%% RECEIVE STRUCTURE
% Specify Receive structure array.
Receive(1).Apod = ones(1,128); % if 64ch Vantage, = [ones(1,64) zeros(1,64)];
Receive(1).startDepth = 0;
Receive(1).endDepth = 200;      % becomes ceil(200*2*4/128)*128/8 = 208
Receive(1).TGC = 1; % Use the first TGC waveform defined above (Time Gain Control)
Receive(1).mode = 0;
Receive(1).bufnum = 1;
Receive(1).framenum = 1;
Receive(1).acqNum = 1;
Receive(1).sampleMode = 'NS200BW';
Receive(1).LowPassCoef = [];
Receive(1).InputFilter = [];


%% PROCESSING OF THE DATA
% Specify an external processing event.
Process(1).classname = 'External';
Process(1).method = 'myProcFunction';
Process(1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                            'srcbufnum',1,...
                            'srcframenum',1,...
                            'dstbuffer','none'};


%% EVENTS                        
% event 1 -- acquisition of RF data
Event(1).info = 'Acquire RF Data.';
Event(1).tx = 1; % use 1st TX structure.
Event(1).rcv = 1; % use 1st Rcv structure.
Event(1).recon = 0; % no reconstruction.
Event(1).process = 0; % no processing
Event(1).seqControl = [1,2];
% sequences for first event
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 50000;
SeqControl(2).command = 'transferToHost';

% event 2 -- processing
Event(2).info = 'Call external Processing function.';
Event(2).tx = 0; % no TX structure.
Event(2).rcv = 0; % no Rcv structure.
Event(2).recon = 0; % no reconstruction.
Event(2).process = 1; % call processing function
Event(2).seqControl = [3,4,5]; % wait for data to be transferred
% sequences event 2
SeqControl(3).command = 'waitForTransferComplete';
SeqControl(3).argument = 2;
SeqControl(4).command = 'markTransferProcessed';
SeqControl(4).argument = 2;
SeqControl(5).command = 'sync';

% event 3 -- jump, do nothing except for jump to event 1
Event(3).info = 'Jump back to Event 1.';
Event(3).tx = 0; % no TX structure.
Event(3).rcv = 0; % no Rcv structure.
Event(3).recon = 0; % no reconstruction.
Event(3).process = 0; % no processing
Event(3).seqControl = 6; % jump back to Event 1
% sequence event 3
SeqControl(6).command = 'jump';
SeqControl(6).argument = 1;


%% UI CONTROL
% Create UI control for channel selection
nr = Resource.Parameters.numRcvChannels;
UI(1).Control = {'UserB1','Style','VsSlider',...
                    'Label','Plot Channel',...
                    'SliderMinMaxVal',[1,128,32],...
                    'SliderStep', [1/nr,8/nr],...
                    'ValueFormat', '%3.0f'};
UI(1).Callback = text2cell('%CB#1');


%% EXTERNAL PLOT FUNCTION
% Create External function for plotting channel data
EF(1).Function = text2cell('%EF#1');


%% SAVE ALL DATA TO MAT FILE FOR VSX
% Save all the structures to a .mat file.
save('L11-4vAcquireRF');

return % Place this return to prevent executing code below


%% DEFINITION OF FUNCTIONS (UI AND EF)

%CB#1
assignin('base', 'myPlotChnl', round(UIValue));
%CB#1


%EF#1
myProcFunction(RData)
persistent myHandle
% If ‘myPlotChnl’ exists, read it for the channel to plot.
if evalin('base','exist(''myPlotChnl'',''var'')')
channel = evalin('base','myPlotChnl');
else
channel = 32; % Channel no. to plot
end
% Create the figure if it doesn’t exist.
if isempty(myHandle)||~ishandle(myHandle)
figure;
myHandle = axes('XLim',[0,1500],'YLim',[-16384 16384], ...
'NextPlot','replacechildren');
end
% Plot the RF data.
plot(myHandle,RData(:,channel));
drawnow
%EF#1