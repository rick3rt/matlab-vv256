% Copyright 2001-2014 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% File name SetUpLDualXdcr_HVMux.m:
% Example script to illustrate the use of two separate transducers imaging
% simultaneously on Vantage-256 system.  For this example, an L12-3v
% transducer is connected to the left connector and a L7-4 to the right
% connector.  This script is a combination of the "L12-3vFlash" and a
% "L7-4Flash" example scripts, with both of them running in parallel but
% independently of each other in the same script.  To accomplish this, we
% define separate structures for SFormat, PData, Recon, and ImageDisplay
% for each transducer.  The script requires a single Trans structure,
% however, defined as if it were a aingle 256 element transducer using both
% connectors- but within that shared structure elements 1:128 represent the
% L12-3v and elements 129:256 represent the L7-4.  The TX and Receive .Apod
% arrays for each transducer will have zeros in the channel positions
% representing the other transducer.  Note however that since both
% transducers are sharing a single Trans structure, they must also share
% the same receive acquisition and processing sample rate as set by
% Trans.frequency.  For this example we are using 4.5 MHz, a reasonable
% compromise for the L12-3v and L7-4. Note also that the GUI controls will generally
% modify settings for the L12-3v only, since it is indexed first in the relevant structures.
% Additional GUI control object(s) can be created to switch between probes.
%
% Two separate TPC profiles are used, to allow independent control of the
% transmit voltage for each transducer using P1 and P2 HV sliders.
%
% Throught this script, comments starting with: % ***Dual Xdcr***
% have been inserted to illustrate the changes made for the dual-imaging
% example

clear all

if exist('Resource', 'var'), Vantage = 1; else Vantage = 0; end

% --- Commonly Changed User Parammeters -------------------------
filename = ('DualXdcr_HVMux'); % name of matfile
numRFFrames = 10; % RF Data frames
frameRateFactor = 5; % Factor for converting sequenceRate to frameRate.
Fc = 4.5; % Trans.frequency
% ---------------------------------------------------------------

% Specify system parameters.
Resource.Parameters.numTransmit = 256; % ***Dual Xdcr*** number of transmit channels.
Resource.Parameters.numRcvChannels = 256; % number of receive channels.
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% ***Dual Xdcr*** We must explicitly specify a connector
% value of zero, to inform the system that both connectors will be used to
% support an overall "transducer" with 256 elements
Resource.Parameters.connector = 0;
% ***Dual Xdcr*** Also set the 'fakeScanhead' parameter, to allow this
% script to run on the HW system with no transducer actually connected and
% avoid confusion when the two distinct transducers are present.
Resource.Parameters.fakeScanhead = 1;

% Specify Trans structure array.
% ***Dual Xdcr*** First we will use computeTrans to define two separate
% Trans1 and Trans2 structures, and then we will use parameters from
% both of those to create the shared Trans structure that will actually be
% used by the script.

% ***Dual Xdcr*** Create the Trans1 structure
Trans1.name = 'L12-3v';
Trans1.units = 'wavelengths';
Trans1 = computeTrans(Trans1); % L12-3v transducer is 'known' transducer so we can use computeTrans.
% ***Dual Xdcr*** Create the Trans2 structure
Trans2.name = 'L7-4';
Trans2.units = 'wavelengths';
Trans2 = computeTrans(Trans2);

% ***Dual Xdcr*** Now use Trans1 and Trans2 to create the shared Trans structure ***
Trans.name = 'custom'; % Must be 'custom' to prevent confusion from the two unique transducer ID's that will actually be connected
Trans.units = 'wavelengths';
Trans.id = hex2dec('0000'); % Dummy ID to be used with the 'fake schanehad' feature
Trans.frequency = Trans1.frequency; % ***Dual Xdcr*** This is the shared frequency to be used by both L12-3v and L7-4
Trans.type = 0; % 1D straight array geometry applies to both L12-3v and L7-4
Trans.numelements = Trans1.numelements + Trans2.numelements; % total over both connectors
% Concatenate the two element position arrays
Trans.ElementPos = [Trans1.ElementPos; Trans2.ElementPos];
% Use a compromise average value for lens correction
Trans.lensCorrection = Trans1.lensCorrection;
% For the following parameters just copy the L12-3v values
Trans.spacing = Trans1.spacing;
Trans.elementWidth = Trans1.elementWidth;
Trans.ElementSens = Trans1.ElementSens;
% For the following use an appropriate shared value
Trans.impedance = 50;
Trans.maxHighVoltage = 50; % set maximum high voltage limit for pulser supply.
Trans.connType = 1;
Trans.HVMux = Trans1.HVMux;
Trans.HVMux.Aperture = [Trans1.HVMux.Aperture; repmat((129:1:256)', 1, size(Trans1.HVMux.Aperture, 2))];
% [~, ~, personalityId] = vdasScanheadQuery();
% Trans.id = personalityId(1);
% Specify TPC structures ... creates two TPC profiles and two HV control sliders.
TPC(1).name = 'L12-3v';
TPC(1).maxHighVoltage = 50;
TPC(2).name = 'L7-4';
TPC(2).maxHighVoltage = 50;

% Specify PData(1) structure array for L12-3v
P(1).startDepth = 5; % Acquisition depth in wavelengths
P(1).endDepth = 192; % This should preferrably be a multiple of 128 samples.

PData(1).PDelta = [Trans1.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P(1).endDepth - P(1).startDepth) / PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans1.numelements * Trans1.spacing) / PData(1).PDelta(1));
PData(1).Size(3) = 1; % single image page
PData(1).Origin = [-Trans1.spacing * (Trans1.numelements - 1) / 2, 0, P(1).startDepth]; % x,y,z of upper lft crnr.

% Specify PData(1) structure array for L7-4
P(2).startDepth = 0; % Acquisition depth in wavelengths
P(2).endDepth = 160; % This should preferrably be a multiple of 128 samples.
PData(2).PDelta = [Trans2.spacing, 0, 0.5];
PData(2).Size(1) = ceil((P(2).endDepth - P(2).startDepth) / PData(2).PDelta(3)); % startDepth, endDepth and pdelta set PData(2).Size.
PData(2).Size(2) = ceil((Trans2.numelements * Trans2.spacing) / PData(2).PDelta(1));
PData(2).Size(3) = 1; % single image page
PData(2).Origin = [-Trans2.spacing * (Trans2.numelements - 1) / 2, 0, P(2).startDepth]; % x,y,z of upper lft crnr.

% Specify Media object.
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2 * 4096; % ***Dual Xdcr*** doubles since there will be two acquisitions per frame, one for each xdcr
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numRFFrames; % e.g., 10 frames used for RF cineloop.
% ***Dual Xdcr***  define the first ImageBuffer and DisplayWindow as usual, for the
% L12-3v
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1; % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
% Resource.ImageBuffer(1).rowsPerFrame = 1024; % this is for maximum depth
% Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L12-3vFlash';
Resource.DisplayWindow(1).pdelta = 0.3;
ScrnSize = get(0, 'ScreenSize');
DwWidth = ceil(PData(1).Size(2) * PData(1).PDelta(1) / Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1) * PData(1).PDelta(3) / Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [50, (ScrnSize(4) - (DwHeight + 150)) / 2, ...% lower left corner position
                                DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1), 0, PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% ***Dual Xdcr*** Now define a second ImageBuffer and DisplayWindow for the
% L7-4
Resource.ImageBuffer(2).datatype = 'double';
Resource.ImageBuffer(2).rowsPerFrame = 1024; % this is for maximum depth
Resource.ImageBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.ImageBuffer(2).numFrames = 10;
Resource.DisplayWindow(2).Title = 'L7-4Flash';
Resource.DisplayWindow(2).pdelta = 0.3;
DwWidth = ceil(PData(2).Size(2) * PData(2).PDelta(1) / Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(2).Size(1) * PData(2).PDelta(3) / Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [600, 200, DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(2).Origin(1), 0, PData(2).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).numFrames = 20;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = gray(256);

% Specify Transmit waveform structure.
% ***Dual Xdcr***
TW(1).type = 'parametric';
TW(1).Parameters = [Trans1.frequency, .67, 2, 1];
TW(2).type = 'parametric';
TW(2).Parameters = [Trans2.frequency, .67, 2, 1];

% ***Dual Xdcr*** Specify TX structure array for the L12-3v.
TX(1).waveform = 1; % use 1st TW structure.
TX(1).aperture = 1;
TX(1).Origin = [0.0, 0.0, 0.0]; % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0, 0.0]; % theta, alpha = 0.
TX(1).Apod = [ones(1, 128), zeros(1, 128)]; % ***Dual Xdcr*** L12-3v uses channels 1:128
TX(1).Delay = computeTXDelays(TX(1));

TX(2) = TX(1);
TX(2).aperture = 65;
TX(2).Delay = computeTXDelays(TX(2));

% ***Dual Xdcr*** Separate TX structure for the L7-4
TX(3).waveform = 2;
TX(3).aperture = TX(2).aperture;
TX(3).Origin = [0, 0, 0]; % set origin to 0,0,0 for flat focus.
TX(3).focus = 0; % set focus to negative for concave TX.Delay profile.
TX(3).Steer = [0, 0];
TX(3).Apod = [zeros(1, 128), ones(1, 128)]; % ***Dual Xdcr*** L7-4 uses channels 129:256
TX(3).Delay = computeTXDelays(TX(3));

% Specify TGC Waveform structure.
TGC(1).CntrlPts = [400, 490, 550, 610, 670, 730, 790, 850];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
TGC(2).CntrlPts = [87, 580, 639, 698, 750, 844, 929, 1023];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
% ***Dual Xdcr*** For our simultaneous dual transducer acquisition scheme,
% we define two interleaved sets of Receive structures.  Each acquisition
% frame will consist of one L12-3vFlash acquisition followed by one L7-4Flash
% acquisition.  This same concept can be easily extended to
% multi-acquisition formats such as FlashAngles, Doppler ensembles, or ray
% line imaging.
maxAcqLength1 = sqrt(P(1).endDepth^2 + (Trans1.numelements * Trans1.spacing)^2) - P(1).startDepth;
maxAcqLength2 = sqrt(P(2).endDepth^2 + (Trans2.numelements * Trans2.spacing)^2) - P(2).startDepth;
wlsPer128 = 128 / (4 * 2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', zeros(1, 256), ...
    'aperture', 1, ...
    'startDepth', P(1).startDepth, ...
    'endDepth', P(1).startDepth + wlsPer128 * ceil(maxAcqLength1 / wlsPer128), ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, 3 * Resource.RcvBuffer(1).numFrames); % ***Dual Xdcr*** two Receive sturctures per frame

% ***Dual Xdcr***  - Set event-specific and transducer-specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % ***Dual Xdcr*** For each frame number i,
    % Receive(3*i-1) and (3*i-2) will be for the L12-3v, and
    % Receive(3*i) will be for the L7-4.
    % Note that most of the L12-3v specific values were used in the initial
    % Receive structure definition so they don't need to be updated here

    Receive(3 * i - 2).Apod(1:96) = 1;
    Receive(3 * i - 2).demodFrequency = TW(1).Parameters(1);
    Receive(3 * i - 2).framenum = i;
    Receive(3 * i - 2).acqNum = 1;
    Receive(3 * i - 2).callMediaFunc = 1;

    Receive(3 * i - 1).Apod(33:128) = 1;
    Receive(3 * i - 1).demodFrequency = TW(1).Parameters(1);
    Receive(3 * i - 1).aperture = 65;
    Receive(3 * i - 1).framenum = i;
    Receive(3 * i - 1).acqNum = 2;

    Receive(3 * i).framenum = i;
    Receive(3 * i).Apod = [zeros(1, 128), ones(1, 128)];
    Receive(3 * i).startDepth = P(2).startDepth;
    Receive(3 * i).endDepth = P(2).startDepth + wlsPer128 * ceil(maxAcqLength2 / wlsPer128);
    Receive(3 * i).TGC = 2;
    Receive(3 * i).acqNum = 3;
    Receive(3 * i).demodFrequency = TW(2).Parameters(1);
    Receive(3 * i).aperture = 65;
end

% Specify Recon structure for L12-3v.
Recon(1) = struct('senscutoff', 0.5, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1, 1], ...
    'ImgBufDest', [1, -1], ...
    'RINums', 1:2);

% Define ReconInfo structure for L12-3v.
ReconInfo = repmat(struct('mode', 'replaceIQ', ...% replace IQ data.
'txnum', 1, ...
    'rcvnum', 1, ...
    'regionnum', 1), 1, 2);
ReconInfo(2).mode = 'accumIQ_replaceIntensity';
ReconInfo(2).txnum = 2;
ReconInfo(2).rcvnum = 2;

% Specify Recon structure for L7-4.
% Copy Recon(1) and then modify values that are different
Recon(2) = Recon(1);
Recon(2).pdatanum = 2;
Recon(2).ImgBufDest = [2, -1];
Recon(2).RINums = 3;

% Define ReconInfo structure for L7-4 using same approach as for Recon.
ReconInfo(3) = ReconInfo(1);
ReconInfo(3).mode = 'replaceIntensity';
ReconInfo(3).txnum = 3;
ReconInfo(3).rcvnum = 3;

% Specify Process structure array for the L12-3v.
pers = 10;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum', 1, ...% number of buffer to process.
'framenum', -1, ...% (-1 => lastFrame)
'pdatanum', 1, ...% number of PData structure to use
'pgain', 1.0, ...% pgain is image processing gain
'reject', 2, ...% reject level
'persistMethod', 'simple', ...
    'persistLevel', pers, ...
    'interpMethod', '4pt', ...
    'grainRemoval', 'none', ...
    'processMethod', 'none', ...
    'averageMethod', 'none', ...
    'compressMethod', 'power', ...
    'compressFactor', 40, ...
    'mappingMethod', 'full', ...
    'display', 1, ...% display image after processing
'displayWindow', 1};

% Specify separate Process structure array for the P6-3.
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum', 2, ...% number of buffer to process.
'framenum', -1, ...% (-1 => lastFrame)
'pdatanum', 2, ...% number of PData structure to use
'pgain', 1.0, ...% pgain is image processing gain
'reject', 2, ...% reject level
'persistMethod', 'simple', ...
    'persistLevel', pers, ...
    'interpMethod', '4pt', ...
    'grainRemoval', 'none', ...
    'processMethod', 'none', ...
    'averageMethod', 'none', ...
    'compressMethod', 'power', ...
    'compressFactor', 40, ...
    'mappingMethod', 'full', ...
    'display', 1, ...% display image after processing
'displayWindow', 2};

% Specify SeqControl structure arrays.

% Set the frame interval using timeToNextAcq
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 50000; % 10000usec = 10msec (~ 100 fps)

% return to matlab for GUI updates every fifth frame
SeqControl(2).command = 'returnToMatlab';

% at the end, jump back and start over
SeqControl(3).command = 'jump';
SeqControl(3).argument = 3; % don't need to repeat the first to events that made initial TPC profile selection

% select TPC profile 2 for the L7-4
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).argument = 2;
SeqControl(4).condition = 'next';

% select TPC profile 1 for the L12-3v
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).argument = 1;
SeqControl(5).condition = 'next';

nsc = 6; % nsc is index of next SeqControl object to be defined

% Specify the Event sequence
n = 1; % n is count of Events

Event(n).info = 'select TPC profile 1 at startup';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc; % set TPC profile command.
n = n + 1;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 1;
SeqControl(nsc).condition = 'immediate';
nsc = nsc + 1;

Event(n).info = 'noop delay for initial profile selection';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc; % noop to allow time for TPC profile transition.
n = n + 1;
SeqControl(nsc).command = 'noop';
SeqControl(nsc).argument = 100e3 / .2; % 100 msec delay
nsc = nsc + 1;

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'acquisition for L12-3v';
    Event(n).tx = 1; % use 1st TX structure.
    Event(n).rcv = 3 * i - 2; % use 1st of the ith pair of Rcv structures.
    Event(n).recon = 0; % no reconstruction.
    Event(n).process = 0; % no processing
    Event(n).seqControl = 1; % time to next acq for frame rate
    n = n + 1;

    Event(n).info = 'acquisition for L12-3v';
    Event(n).tx = 2; % use 1st TX structure.
    Event(n).rcv = 3 * i - 1; % use 1st of the ith pair of Rcv structures.
    Event(n).recon = 0; % no reconstruction.
    Event(n).process = 0; % no processing
    Event(n).seqControl = [4, 1]; % time to next acq for frame rate
    n = n + 1;

    Event(n).info = 'acquisition for L7-4';
    Event(n).tx = 3; % use 1st TX structure.
    Event(n).rcv = 3 * i; % use 2nd of the ith pair of Rcv structures.
    Event(n).recon = 0; % no reconstruction.
    Event(n).process = 0; % no processing
    Event(n).seqControl = [5, 1, nsc]; % use SeqControl struct defined below.
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n + 1;

    Event(n).info = 'Reconstruct & display L12-3v';
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 1; % separate reconstruction for each transducer
    Event(n).process = 1; % separate image display processing for each transducer
    Event(n).seqControl = 0;
    n = n + 1;

    Event(n).info = 'Reconstruct & display L7-4';
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 2; % separate reconstruction for each transducer
    Event(n).process = 2; % separate image display processing for each transducer

    if floor(i / frameRateFactor) == i / frameRateFactor% Exit to Matlab every 5th frame
        Event(n).seqControl = 2; % return to Matlab
    else
        Event(n).seqControl = 0;
    end

    n = n + 1;
end

Event(n).info = 'Jump back to third event to repeat the sequence';
Event(n).tx = 0; % no TX
Event(n).rcv = 0; % no Rcv
Event(n).recon = 0; % no Recon
Event(n).process = 0;
Event(n).seqControl = 3; % jump command

% Save all the structures to a .mat file.
save(filename);
%  VSX % run automatically

return
