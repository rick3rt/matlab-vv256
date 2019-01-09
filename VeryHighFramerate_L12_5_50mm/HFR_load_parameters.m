function P = HFR_load_parameters(acqMode, simMode)
% HFR_LOAD_PARAMETERS: Function that load the global parameters that can be
% used for Acquisition Modes:
%   'scan'    low frame rate scan
%   'acq'     high frame rate imaging
%   'recon'   reconstruction of IQ data
%
% Simulate Modes:
%   0       actual imaging with system
%   1       simulation mode
%
% Example:
%   HFR_LOAD_PARAMETERS('scan', 0)   
%       loads parameters for high frame rate
%       imaging with the system
%

%% User variables

% set number of frames
numOfFrames_scan = 30;          % scanning is low frame rate 'exploring'

numOfFrames_imaging = 1000;      % imaging is high frame rate acq

% set imaging speed
fs = 100; % Hz - imaging frequency


% imaging depth
startDepth = 2; % Acquisition depth in wavelengths
endDepth = 128; % This should preferrably be a multiple of 128 samples.



%% check if modes are arguments.
if ~exist('acqMode','var')
    acqMode = 'scan';
end

if ~exist('simMode','var')
    simMode = 1; % default to sim mode, not real imaging
end


%% constant for all modes
P.acqMode = acqMode;

P.startDepth = startDepth; % Acquisition depth in wavelengths
P.endDepth = endDepth; % This should preferrably be a multiple of 128 samples.

P.imageRate = fs;

P.totalFrames = numOfFrames_imaging;
P.numApertures = 3;
P.numFramesPerSuperFrame = 30; % limited by memory in Verasonics Vantage (?) LOOKUP
P.numSuperFrames = ceil(numOfFrames_imaging/P.numFramesPerSuperFrame);
P.numAcqs = P.numApertures*P.numFramesPerSuperFrame;    % no. of Acquisitions in a Receive frame (this is a "superframe")
P.totalFrames = P.numFramesPerSuperFrame*P.numSuperFrames;

% sim mode
P.simMode = 1;      % set to acquire data using Vantage 128/256 hardware

% timing, duration of timeToNextAcq
% can range from 10 – 4190000 microsec
P.time.aperture = 200;              % usec
P.time.frame = 1/fs*1e6;            % usec
P.time.superFrame = P.time.frame;   % usec
P.time.unit = 'microseconds';

% Image processing
P.persistLevel = 20;

% not sure what they do... but same for all modes
P.frameRateFactor = 2; 
P.persistLevel = 20; % for image display (Processing function)


%% mode dependent
switch acqMode
    
    case 'scan'
        % LFR - scanning
        P.numFramesPerSuperFrame = 30; 
        P.numSuperFrames = 1;
        P.numAcqs = P.numApertures*P.numFramesPerSuperFrame;
        P.totalFrames = P.numFramesPerSuperFrame*P.numSuperFrames;
        
        P.simMode = simMode;
        
        fs = 25;
        P.time.aperture = 200;              % usec
        P.time.frame = 1/fs*1e6;            % usec
        P.time.superFrame = P.time.frame;   % usec
        
    case 'acq'
        % HFR - imaging, acquiring RF data
        P.simMode = simMode;

    case 'recon'
        % IQR - reconstructing IQData
        P.simMode = 2;
        
    otherwise
        % none of the above, i.e. invalid
        warning('Unexpected imaging mode. No parameters returned.')
        P = struct();
end


end