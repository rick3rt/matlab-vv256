% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name: SetUpL22_14vFlashAnglesInterleave.m 
%               - Example of plane wave imaging at 83.3 MHz with "interleaved sampling" using two 41.7 MHz acquisitions
%
% Description: 
%   Sequence programming file for L22-14v Linear array, using plane wave
%   transmits with multiple steering angles using 4X sampling. All 128 
%   transmit and receive channels are active for each acquisition. 
%   This script demonstrates the "interleaved sampling" mechanism for RF 
%   data acquisition.  In this script, the L22-14v is operated at a nominal
%   center frequency for Recon processing of 20.83 MHz, with 4X sampling so
%   the required RF data acquisition sample rate is 83.3 MHz. Since the HW
%   system cannot acquire data at sample rates above 62.5 MHz, we acquire 
%   two sets of receive data for each line, at an A/D sample rate of 41.7 MHz.
%   These pairs of acquisition data lines are then interleaved in an external
%   processing function, to produce a single composite line sampled at 83.3
%   MHz as expected by Recon.  The data acquired in the second acquistion of
%   each pair must have its data shifted by half a sample period, so it can
%   be interleaved with the first line.  This is accomplished by adding an
%   offset to the transmit delays for the second line equal to that
%   half-sample interval; this in effect shifts each receive sample's data
%   earlier by the desired one half sample period(other than the interleaved
%   acquisition at 83.3 MHz, this script is identical to L22-14vFlashAngles
%   which uses acquisition at 62.5 MHz with no interleave.
%
% Last update:
% 12/07/2015 - modified for SW 3.0

clear all
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

RcvProfile.LnaZinSel = 31;

na = 7;      % Set na = number of angles.
if (na > 1)
    dtheta = (36*pi/180)/(na-1); 
    P.startAngle = -36*pi/180/2; 
else
    dtheta = 0; 
    P.startAngle=0; 
end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 25; % mfr data sheet lists 30 Volt limit

% Specify PData structure array.
PData(1).PDelta = [0.4, 0, 0.25];
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
Resource.RcvBuffer(1).rowsPerFrame = na*4096*4; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L22-14v FlashAngles Interleave, 4X sampling at 83 MHz';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.  
TW.type = 'parametric';
TW.Parameters = [18, 0.67, 3, 1];   % 18 MHz center frequency, 67% pulsewidth 1.5 cycle burst

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2*na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end

% TX(1) delay will be 1/4 of the wls, based on the demodulation frequency
% TX.delay will be converted to time based on the Trans.frequency
for n = 1:2:2*na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)/2*dtheta),0.0];
    TX(n+1).Steer = TX(n).Steer;
    TX(n+1).Delay = computeTXDelays(TX(n+1));
    TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/20.83); 
end

% with 4X sampling, the 20.83 MHz transmit frequency maps to an acquisition
% sample rate of 83.3 MHz, which is 250 MHz divided by 3.

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,511,716,920,1023,1023,1023,1023]; %[400,550,650,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays. 
% - We need 2*na Receives for every acquisition frame, since two
% acquisitions are needed to produce one interleaved acquisition line of
% data.  An addition na "Dummy" Receive structures are needed to convey to
% Recon the actual interleaved line length.

% For the interleaved acquisition approach with 2X interleave, we must take
% into account the actual sample rate at which the digital filters in the
% CGD FPGA will be operating: twice the center frequency set by
% Trans.frequency, not the typical 4X factor. This means that the Nyquist
% limit for the CGD filters will be at Trans.frequency; the higher half of
% the transducer signal frequency spectrum will be folded over the lower
% half due to aliasing.  (The subsequent interleave combination of the two
% acquisition events will unfold this aliasing).  Therefore the filter
% actually used for the CGD Input Filter must be defined as a high-pass
% filter;  the net effect after interleave will be a symmetric bandpass
% filter centered at 20.83 MHz.

% Four highpass coefficient arrays have been defined below, representing
% fractional bandwidths of 100%, 67%, 50%, and 25% relative to Fc at 20.83
% MHz.  The coefficient arrays listed below were developed using the matlab
% fdatool, with a Hamming window except for "HighPassCoef100K" using a
% Kaiser window.

% Note that for the L22-14v we would ideally use a bandpass filter centered
% near the transducer's nominal center frequency of 18 MHz, not the 20.83
% MHz forced on the CGD filters.  In the external processing function that
% does the interleave to produce the 83.3 MHz sampled RF data, a filter
% could be added operating at that rate to allow any desired center
% frequency.
HighPassCoef100K =[-0.0000    0.0161   -0.0000   -0.0182   -0.0000    0.0208   -0.0000...
                   -0.0242   -0.0000    0.0288   -0.0000   -0.0355   -0.0000    0.0458...
                   -0.0000   -0.0644   -0.0000    0.1076   -0.0000   -0.3231    0.5076];
               
HighPassCoef100 = [-0.0000    0.0014   -0.0000   -0.0024   -0.0000    0.0046   -0.0000...
                   -0.0081   -0.0000    0.0136   -0.0000   -0.0217   -0.0000    0.0341...
                   -0.0000   -0.0551   -0.0000    0.1009   -0.0000   -0.3169    0.5006];
                
HighPassCoef67 = [  0.0011   -0.0012   -0.0000    0.0021   -0.0029   -0.0000    0.0053...
                   -0.0070   -0.0000    0.0117   -0.0148   -0.0000    0.0235   -0.0294...
                   -0.0000    0.0476   -0.0627   -0.0000    0.1345   -0.2735    0.3326];                
                
HighPassCoef50 = [ -0.0000   -0.0010    0.0018   -0.0017   -0.0000    0.0032   -0.0061...
                    0.0057   -0.0000   -0.0096    0.0171   -0.0153   -0.0000    0.0240...
                   -0.0429    0.0389   -0.0000   -0.0711    0.1552   -0.2233    0.2494];

HighPassCoef25 = [  0.0013   -0.0013    0.0013   -0.0009   -0.0000    0.0017   -0.0043...
                    0.0075   -0.0105    0.0125   -0.0121    0.0083   -0.0000   -0.0130...
                    0.0303   -0.0508    0.0724   -0.0929    0.1097   -0.1208    0.1247];
                
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
% ensure half sample rate lines are an integer number of blocks, and then
% double for full sample rate lines.
linelength = 2*maxAcqLength;
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', linelength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'demodFrequency', 20.83,... % allowed frequency 
                        'InputFilter', HighPassCoef100, ... 
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 2*na*Resource.RcvBuffer(1).numFrames);
                    
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(2*na*(i-1)+1).callMediaFunc = 1;
    for j = 1:na
        % pairs of acquisition receives to buffer 1 for interleaved
        % acquisition at 1/2 sample rate:
        % First of each pair
        Receive(2*na*(i-1)+2*j-1).framenum = i;
        Receive(2*na*(i-1)+2*j-1).acqNum = 2*j-1;
        % second of each pair
        Receive(2*na*(i-1)+2*j).framenum = i;
        Receive(2*na*(i-1)+2*j).acqNum = 2*j;        
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame. 
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na  
        ReconInfo(j).txnum = 2*j-1; % only the first TX is required for ReconInfo
        ReconInfo(j).rcvnum = 2*j-1; 
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
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

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - (na-1)*160;  % 20 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                      % Acquire frame
        Event(n).info = 'Acquire angle, 1st half interleave samples';
        Event(n).tx = 2*j-1;   
        Event(n).rcv = 2*na*(i-1)+2*j-1;   
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = 2;
        n = n+1;
        
        Event(n).info = 'Acquire angle, 2nd half interleave samples';
        Event(n).tx = 2*j;   
        Event(n).rcv = 2*na*(i-1)+2*j;   
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = 2;
        n = n+1;
        
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 1;      
    Event(n).process = 1;    
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame 
        Event(n).seqControl = 4; 
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
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
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L22-14vFlashAnglesInterleave');
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
