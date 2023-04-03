function [Fastro, Fback, Trel] = GetAstroFluor( expt, runInfo, projParam, backROI, cellROI, somaROI, procROI, varargin )

IP = inputParser;
addRequired( IP, 'expt', @isstruct) %dataDir = 'D:\2photon\AB18\190905\008\';
addRequired( IP, 'runInfo', @isstruct) 
addRequired( IP, 'projParam', @isstruct)
addRequired( IP, 'backROI', @isstruct)
addRequired( IP, 'cellROI', @isstruct)
addRequired( IP, 'somaROI', @isstruct)
addRequired( IP, 'procROI', @isstruct)
addParameter( IP, 'refTime', [], @isdatetime)
%addParameter( IP, 'save', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'show', false, @islogical )
parse( IP, expt, runInfo, projParam, backROI, cellROI, somaROI, procROI, varargin{:} );
refTime = IP.Results.refTime;
%setSave = IP.Results.save;
overwrite = IP.Results.overwrite;
show = IP.Results.show;
Ncell = numel(cellROI);
savePath = sprintf('%s%s_astroFluor.mat', expt.dir, expt.name);
tic
if ~exist(savePath,'file') || overwrite
    if isempty( refTime ), refTime = runInfo(1).timestamp; end
    % Load concatenated green projection
    greenCatPath = projParam.path.cat.reg.z{2};
    fprintf('\nLoading %s', greenCatPath)
    greenCat = loadtiff(greenCatPath);
    % EXTRACT ROI FLUORESCENCE FROM MOVIE
    % get pixel-level fluorescence for each ROI
    fprintf('\nExtracting fluorescence:  '); 
    FcellPix = cell(1,Ncell); FsomaPix = cell(1,Ncell); FprocPix = cell(1,Ncell);
    for c = flip(1:Ncell)
        for p = flip(1:cellROI(c).Npix), FcellPix{c}(:,p) = greenCat( cellROI(c).pix(p,2),  cellROI(c).pix(p,1), : ); end
        for p = flip(1:somaROI(c).Npix), FsomaPix{c}(:,p) = greenCat( somaROI(c).pix(p,2),  somaROI(c).pix(p,1), : ); end
        for p = flip(1:procROI(c).Npix), FprocPix{c}(:,p) = greenCat( procROI(c).pix(p,2),  procROI(c).pix(p,1), : ); end
    end
    for p = flip(1:backROI.Npix), FbackPix(:,p) = greenCat( backROI.pix(p,2),  backROI.pix(p,1), : ); end %#ok<AGROW>
    % calculate mean fluorescence
    Trel = cell(1,expt.Nruns);  Fback = cell(1,expt.Nruns);
    Fastro = repmat( struct('cell',[], 'soma',[], 'proc',[]), 1, expt.Nruns );
    for runs = 1:expt.Nruns
        Trel{runs} = (projParam.Tproj{runs}-projParam.Tproj{runs}(1) + seconds(runInfo(runs).timestamp - refTime))/60; % shift time to be relative to reference time, convert to minutes
        runBins = projParam.binLims(runs)+1:projParam.binLims(runs+1);
        Fback{runs} = mean( FbackPix(runBins,:), 2, 'omitnan');
        Fback_med(runs) = median(Fback{runs});
        for c = flip(1:Ncell)
            Fastro(runs).cell(:,c) = mean( FcellPix{c}(runBins,:), 2, 'omitnan');
            Fastro(runs).soma(:,c) = mean( FsomaPix{c}(runBins,:), 2, 'omitnan');
            Fastro(runs).proc(:,c) = mean( FprocPix{c}(runBins,:), 2, 'omitnan');
        end
        % Background subtraction
        Fastro(runs).cell_sub = Fastro(runs).cell - Fback{runs} + Fback_med(runs);
        Fastro(runs).soma_sub = Fastro(runs).soma - Fback{runs} + Fback_med(runs);
        Fastro(runs).proc_sub = Fastro(runs).proc - Fback{runs} + Fback_med(runs);

        %{
        figure
        for c = 1:Ncell
            subplot(2,1,1); plot(Fastro(runs).cell);
            subplot(2,1,2); plot(Fastro(runs).cell_sub);
        end
        %}

        % Normalization
        % Calculate baseline signals Fo using rolling percentile
        %{
        fprintf('\nCalculating baseline signal');
        Fo = nan(size(F));
        for run = 1:Nruns
            Fastro(runs).cell_base = MovingPercentile(Fastro(runs).cell_sub, 10, windowSize, 'pre');
            Fo(runScans{run},:) = MovingPercentile(F(runScans{run},:), basePrct, windowSize, 'pre'); % high-pass filter movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
        end
        Fotot = MovingPercentile(Ftot, basePrct, windowSize, 'pre'); 
        dFF = (F-Fo)./Fo; % (Ffilt-Fo)./Fo; %dFF(tempCC.PixelIdxList{c},:) = (F-Fo)./Fo;
        dFFtot = (Ftot-Fotot)./Fotot;

        %F = Froi_lp - Fnp_lp + mean(Fnp_lp,1,'omitnan');  % Froi - Fnp + mean(Fnp,1,'omitnan'); % plot( F(:,1) );
        %Ftot = Fall_lp - Fback_lp + mean(Fback_lp,1,'omitnan');
        %}
    end
    % Save the results
    fprintf('\nSaving %s', savePath);
    save(savePath, 'expt', 'projParam', 'refTime', 'Ncell','Fastro', 'Trel', 'Fback', 'Fback_med' );
    toc
else
    fprintf('\nLoading %s', savePath);
    load(savePath, 'Fastro', 'Trel', 'Fback', 'Ncell');
    toc
end

% Plot the results (optional)
if show
    cellColor = distinguishable_colors(Ncell);
    figure('WindowState','max')
    for runs = 1:expt.Nruns
        %plot( Trel{runs}, Fback{runs}, 'k'); hold on;
        for c = 1:Ncell
            plot( Trel{runs}, Fastro(runs).cell_sub(:,c), 'color', cellColor(c,:)); hold on;
        end
    end
    %{
    h(1) = plot( T, Fback, 'k' ); hold on;
    h(2) = plot( T, Fastro.cell, 'b' );
    h(3) = plot( T, Fastro.soma, 'r' );
    h(4) = plot( T, Fastro.proc, 'g' );
    legend(h, {'Background','Cell','Soma','Processes'});
    %}
    xlim([-Inf,Inf]);
    xlabel('Time (min)'); ylabel('Raw Fluorescence');
end
end