function [speed, instSpeed, accel, binSpeed, speedBout, speedFrame, speedClean, longStill] = MakeSpeed( metadata, T, quadData, varargin )
%MakeSpeed Return the speed of a mouse running on a 3d printed wheel from date from the quadrature encoder
    
IP = inputParser;
addRequired( IP, 'metadata', @isstruct )
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'quadData', @isnumeric )
addParameter( IP, 'gaussWidth', 0.5, @isnumeric )
addParameter( IP, 'gaussSigma', 0.1, @isnumeric )
addParameter( IP, 'minDur', 0.3, @isnumeric ) % minimum duration, in seconds, for a bout to not be immediately deleted
addParameter( IP, 'minPeak', 3, @isnumeric ) % minimum peak speed for a bout to not be immediately deleted
addParameter( IP, 'maxGap', 0.5, @isnumeric ) % merge 2 bouts which end/begin within this many seconds of eachother
addParameter( IP, 'finDur', 1, @isnumeric ) % minimum length, in seconds, for FINAL bouts to not be deleted
addParameter( IP, 'show', false, @islogical )
parse( IP, metadata, T, quadData, varargin{:} ); 
gaussWidth = IP.Results.gaussWidth; 
gaussSigma = IP.Results.gaussSigma; 
minDur = IP.Results.minDur;
minPeak = IP.Results.minPeak;
maxGapDur = IP.Results.maxGap;
finDur = IP.Results.finDur; 
show = IP.Results.show; 

wheel_diameter = 14; % in cm
wheel_tabs = 44; 
wheel_circumference = wheel_diameter*pi;
step_size = wheel_circumference/(wheel_tabs*2);

instSpeed = zeros(length(quadData), 1);
if ~isempty(instSpeed)
    instSpeed(2:end) = diff(quadData);
    instSpeed(2) = 0;
    instSpeed = instSpeed*step_size*metadata.rate; % fps
end
instSpeed = double(instSpeed);

% Filter the inst speed
gaussFilt = MakeGaussFilt( gaussWidth, 0, gaussSigma, metadata.rate );
speed = filtfilt( gaussFilt, 1, instSpeed ); %conv( instantaneous_speed, gaussFilt, 'same' ); % conv(instantaneous_speed, ones(ceil(framerate/4), 1)/ceil(framerate/4), 'same');


% Plot
%{
if show
    yyaxis right
    plot( instantaneous_speed ); hold on;
    ylabel('Instantaneous Speed (cm/s)');
    yyaxis left
    plot( speed, 'k', 'LineWidth',2 );
    xlabel('Frame'); ylabel('Averaged Speed (cm/s)');
end
%}
if nargout > 2
    % Calculate acceleration
    accel = nan(length(quadData),1); 
    dT = 1/metadata.rate; %median(diff(T)); %
    accel(2:end) = diff(speed)/dT;

    % Detect sustained bouts of quadData
    [binSpeed, speedBout, ~] = speedBouts(T, speed, accel, 'minFrame', round(minDur*metadata.rate), 'minPeak', minPeak, 'mergeSep', round(maxGapDur*metadata.rate), ...
        'finMinFrm', round(finDur*metadata.rate), 'show',show, 'pause',false ); 
    speedFrame.run = find(binSpeed); %[speedBout.frames];  
    speedFrame.still = find(~binSpeed); %1:Nframe; speedFrame.still(speedFrame.run) = [];

    speedClean = speed; speedClean(speedFrame.still) = 0;

    % Find the longest periods of stillness
    tempConn = bwconncomp( ~binSpeed );
    stillConn = tempConn.PixelIdxList;
    [~, sortInd] = sort( cellfun(@length, stillConn), 'descend');
    stillConn = stillConn( sortInd  );
    longStill.frame = stillConn{1};
    longStill.dur = diff(T(longStill.frame([1,end])));
end
end