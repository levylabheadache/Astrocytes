function [Tmid, binSumm, Nbin, binInt] = BinCombinedData(Tcomb, fovComb, Fcomb, dFFcomb, Acomb, locoComb, varargin)
% Examine the combined data at uniform time bins around the reference time
IP = inputParser;
addRequired( IP, 'Tcomb', @isnumeric )
addRequired( IP, 'fovComb', @isnumeric )
addRequired( IP, 'Fcomb', @isstruct )
addRequired( IP, 'dFFcomb', @isstruct )
addRequired( IP, 'Acomb', @isstruct )
addRequired( IP, 'locoComb', @isnumeric )
%addOptional( IP, 'cellEvt', {}, @iscell )
addParameter( IP, 'Trange', [-60,120], @isnumeric ) % how many minutes before/after reference time should we bin?
addParameter( IP, 'interval', 1, @isnumeric ) % length of bin interval (minutes)
addParameter( IP, 'minFrac', 0.75, @isnumeric ) % for a bin to be included, must have observed at least minFrac percent of the total bin duration
addParameter( IP, 'frameSec', 3/15.49, @isnumeric ) % seconds per frame, factor of 3 due to frame-binning for denoising
addParameter( IP, 'show', false, @islogical ) % true
parse( IP, Tcomb, fovComb, Fcomb, dFFcomb, Acomb, locoComb, varargin{:} ); % , velocity
Trange = IP.Results.Trange;
binInt = IP.Results.interval;
minFrac = IP.Results.minFrac;
frameSec = IP.Results.frameSec;
show = IP.Results.show;
if any(~isfinite(Trange)) || Trange(1) >= Trange(2), error('Trange improperly set'); end
% Set binning variables and parameters
Tbin = Trange(1):binInt:Trange(2);
Tmid = Tbin(1:end-1) + binInt/2; 
Nbin = numel(Tmid);
minObs = minFrac*binInt;
Nastro = size(Fcomb.cell,2);
blankStruct = struct('cell',nan(1,Nastro), 'soma',nan(1,Nastro), 'proc',nan(1,Nastro));
binStruct = struct('fov',NaN, 'Nframe',NaN, 'obs',NaN, 'good',false, 'Nastro',NaN, 'topSpeed',NaN, 'travel',NaN, 'Nbout',NaN, 'runTime',NaN, ...
    'Fmed',blankStruct,  'Fcen',blankStruct, 'Abin',blankStruct, 'dFFmed',blankStruct, 'SPratio',nan(1,Nastro), 'LFratio',nan(1,Nastro));
binSumm = repmat( binStruct, Nbin, 1);
% Set Chronux parameters
chronuxParam.Fs = 1/frameSec;
chronuxParam.err = 0;
chronuxParam.pad = 0;
lowFreq = [0, 0.2];
freqRes = 0.1; % Hz 

tic;
for t = 1:Nbin
    binFrames = find( Tcomb(:,1) >= Tbin(t) & Tcomb(:,1) < Tbin(t+1) );
    % Which FoV are observed during this bin? If more than one, use only the most-observed FoV.
    tempFOV = fovComb(binFrames);
    binFov  = unique(tempFOV)';  Nfov = numel(binFov);
    if Nfov > 1
        %fprintf('\nt = %i: multiple FoV observed', t );
        fovFrame = cell(1,Nfov);
        for f = 1:Nfov, fovFrame{f} = find( tempFOV == binFov(f) ); end
        NfovFrame = cellfun( @numel, fovFrame );
        [~,fMax] = max( NfovFrame );
        binFrames = binFrames( fovFrame{fMax} );
        binSumm(t).fov = binFov(fMax); 
    else
        binSumm(t).fov = binFov;
    end
    binSumm(t).Nframe = numel(binFrames);
    if binSumm(t).Nframe > 0
        binSumm(t).obs = Tcomb(binFrames(end),1) - Tcomb(binFrames(1),1);  % dT = 
        if binSumm(t).obs >= minObs
            binSumm(t).good = true;
            % Locomotive properties
            binSumm(t).topSpeed = max( abs(locoComb(binFrames,1)) ); %
            binSumm(t).travel = frameSec*sum(abs(locoComb(binFrames,1)));          
            binSumm(t).runTime = frameSec*sum(abs(locoComb(binFrames,1)) > 0);  %numel(find(locoComb(binFrames,1))); 
            if any(locoComb(binFrames,2))
                binSumm(t).Nbout = numel( unique(locoComb(binFrames,2) > 0) ); 
            end
            % Unnormalized fluor properties
            binSumm(t).Fmed.cell = median( Fcomb.cell(binFrames,:), 1 );
            binSumm(t).Fmed.soma = median( Fcomb.soma(binFrames,:), 1 );
            binSumm(t).Fmed.proc = median( Fcomb.proc(binFrames,:), 1 );
            binSumm(t).SPratio = binSumm(t).Fmed.soma/binSumm(t).Fmed.proc;
            binSumm(t).Nastro = sum(~isnan(binSumm(t).Fmed.cell), 2); % How many astrocytes were observed in this bin?
            % Run-censored median fluorescence
            censBinFrame = binFrames( locoComb(binFrames,3) ~= 1 ); % Find periRun frames for censoring
            if numel(censBinFrame)*frameSec > 5
                binSumm(t).Fcen.cell = median( Fcomb.cell(censBinFrame,:), 1 );
                binSumm(t).Fcen.soma = median( Fcomb.soma(censBinFrame,:), 1 );
                binSumm(t).Fcen.proc = median( Fcomb.proc(censBinFrame,:), 1 );
            end
            % Normalized fluor properties
            binSumm(t).dFFmed.cell = median( dFFcomb.cell(binFrames,:), 1 );
            binSumm(t).dFFmed.soma = median( dFFcomb.soma(binFrames,:), 1 ); 
            binSumm(t).dFFmed.proc = median( dFFcomb.proc(binFrames,:), 1 );
            
            % Activity per time
            binSumm(t).Abin.cell = sum( Acomb.cell(binFrames,:), 1 )/binSumm(t).obs;
            binSumm(t).Abin.soma = sum( Acomb.soma(binFrames,:), 1 )/binSumm(t).obs; 
            binSumm(t).Abin.proc = sum( Acomb.proc(binFrames,:), 1 )/binSumm(t).obs;
            
            % Low-frequency fraction
            timeRes = binSumm(t).obs*60; % convert from minutes to seconds
            TW = timeRes*freqRes;
            K = floor(2*TW-1); 
            chronuxParam.tapers = [TW, K]; 
            binSumm(t).LFratio = LowFreqRatio(dFFcomb.cell(binFrames,:), chronuxParam, lowFreq);               
        end
    end
end
if show
    figure('WindowState','max');
    sp(1) = subplot(3,1,1); 
    plot( Tcomb, Fcomb.cell ); hold on
    Fmed = [binSumm.Fmed];
    FmedCell = vertcat( Fmed.cell ); 
    plot( Tmid, FmedCell, 'LineWidth', 2 ); 
    ylabel('Median Fluorescence'); 
    
    sp(2) = subplot(3,1,2); 
    plot( Tcomb, locoComb(:,1) ); hold on;
    plot( Tmid, [binSumm.topSpeed] );
    legend('Velocity','Bin Top Speed');
    ylabel('Velocity (cm/s)');  
    
    sp(3) = subplot(3,1,3); 
    LFratio = vertcat(binSumm.LFratio);
    plot( Tmid, LFratio ); hold on;
    ylabel('Normalized Low Frequency Power');  
    
    xlim(Tcomb([1,end])); % xlim([-Inf,Inf]);
    xlabel('Time Post-Injection (min)'); 
    linkaxes(sp,'x'); 
end

