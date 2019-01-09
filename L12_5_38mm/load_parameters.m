function Params = load_parameters(acqMode, simMode)
% LOAD_PARAMETERS: Function that load the global parameters that can be used
% for Acquisition Modes:
%   LFR     low frame rate scan
%   HFR     high frame rate imaging
%   IQR     reconstruction of IQ data
%
% Simulate Modes:
%   0       actual imaging with system
%   1       simulation mode
%
% Example:
%   LOAD_PARAMETERS('HFR', 0)   
%       loads parameters for high frame rate
%       imaging with the system
%
%

%% check if modes are arguments.

if ~exist('acqMode','var')
    acqMode = 'LFR';
end

if ~exist('simMode','var')
    simMode = 1; % default to sim mode, not real imaging
end

%% set number of frames
numOfFrames_scan = 30;          % scanning is low frame rate 'exploring'
numOfFrames_imaging = 100;      % imaging is high frame rate recording


%% constant for all modes
Params.acqMode = acqMode;
Params.startDepth = 2;
Params.endDepth = 128;


% not sure what they do... but same for all modes
Params.frameRateFactor = 2; 
Params.persistLevel = 20; % for image display (Processing function)


%% mode dependent
switch acqMode
    
    case 'LFR'
        % LFR - scanning
        Params.numOfFrames = numOfFrames_scan; 
        Params.simMode = simMode;
        
    case 'HFR'
        % HFR - imaging
        Params.numOfFrames = numOfFrames_imaging;
        Params.simMode = simMode;

    case 'IQR'
        % IQR - reconstructing
        Params.numOfFrames = numOfFrames_imaging; 
        Params.simMode = 2;
        
    otherwise
        % none of the above, i.e. invalid
        warning('Unexpected imaging mode. No parameters returned.')
        Params = struct();
end


end