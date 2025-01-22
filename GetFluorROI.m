function [Froi, Fback, Trel] = GetFluorROI( expt, runInfo, projParam, backROI, cellROI, varargin ) % , somaROI, procROI
checkInput = @(x)(isstruct(x) || iscell(x));
IP = inputParser;
addRequired( IP, 'expt', @isstruct) %dataDir = 'D:\2photon\AB18\190905\008\';
addRequired( IP, 'runInfo', @isstruct)
addRequired( IP, 'projParam', @isstruct)
addRequired( IP, 'backROI', checkInput)
addRequired( IP, 'cellROI', checkInput)
addOptional( IP, 'somaROI', [], checkInput)
addOptional( IP, 'procROI', checkInput)
addParameter( IP, 'refTime', [], @isdatetime)
%addParameter( IP, 'save', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'show', false, @islogical )
parse( IP, expt, runInfo, projParam, backROI, cellROI, varargin{:} ); % , somaROI, procROI
somaROI = IP.Results.somaROI;
procROI = IP.Results.procROI;
refTime = IP.Results.refTime;
%setSave = IP.Results.save;
overwrite = IP.Results.overwrite;
show = IP.Results.show;
Nz = numel(cellROI);
Ncell = cellfun(@numel, cellROI); %numel(cellROI);
somaToggle = logical(numel(somaROI));
savePath = sprintf('%s%s_FluorROI.mat', expt.dir, expt.name);
tic
if ~exist(savePath,'file') || overwrite
    if isempty( refTime ), refTime = runInfo(1).timestamp; end
    Fcell_pix = cell(1,Nz); Fsoma_pix = cell(1,Nz); Fproc_pix = cell(1,Nz); Fback_pix = cell(1,Nz);
    Trel = cell(1,expt.Nruns);  Fback = cell(expt.Nruns, Nz); Fback_med = nan(expt.Nruns, Nz);
    Froi = repmat( struct('cell',[], 'soma',[], 'proc',[]), expt.Nruns, Nz);
    for z = 1:Nz
        % Load concatenated green-channel projection
        if expt.Nruns > 1 && ~isempty(projParam.path.cat.reg.z{2,1,z}) %{1,2,z}
            greenCatPath = projParam.path.cat.reg.z{2,1,z}; %{1,2,z}
        elseif expt.Nruns == 1 && ~isempty(projParam.path.run.raw.z{1,2,z}) %projParam.path.run.reg.z{1,2,z})
            greenCatPath = projParam.path.run.raw.z{1,2,z};
        else
            error('Registered projection does not exist!')
        end
        fprintf('\nLoading %s', greenCatPath)
        greenCat = loadtiff(greenCatPath);
        
        % EXTRACT ROI FLUORESCENCE FROM MOVIE
        % get pixel-level fluorescence for each ROI and subROI
        fprintf('\nExtracting fluorescence:  ');
        Fcell_pix{z} = cell(1,Ncell(z)); Fsoma_pix{z} = cell(1,Ncell(z)); Fproc_pix{z} = cell(1,Ncell(z));
        for c = flip(1:Ncell(z))
            for p = flip(1:cellROI{z}(c).Npix), Fcell_pix{z}{c}(:,p) = greenCat( cellROI{z}(c).pix(p,2),  cellROI{z}(c).pix(p,1), : ); end
            if somaToggle
                for p = flip(1:somaROI{z}(c).Npix), Fsoma_pix{z}{c}(:,p) = greenCat( somaROI{z}(c).pix(p,2),  somaROI{z}(c).pix(p,1), : ); end
                for p = flip(1:procROI{z}(c).Npix), Fproc_pix{z}{c}(:,p) = greenCat( procROI{z}(c).pix(p,2),  procROI{z}(c).pix(p,1), : ); end
            end
        end
        for p = flip(1:backROI{z}.Npix), Fback_pix{z}(:,p) = greenCat( backROI{z}.pix(p,2),  backROI{z}.pix(p,1), : ); end
        
        % calculate mean fluorescence       
        for runs = 1:expt.Nruns
            Trel{runs} = (projParam.Tproj{runs}-projParam.Tproj{runs}(1) + seconds(runInfo(runs).timestamp - refTime))/60; % shift time to be relative to reference time, convert to minutes
            runBins = projParam.binLims(runs)+1:projParam.binLims(runs+1);
            Fback{runs,z} = mean( Fback_pix{z}(runBins,:), 2, 'omitnan');
            Fback_med(runs,z) = median(Fback{runs,z}, 'omitnan');
            for c = flip(1:Ncell(z))
                Froi(runs,z).cell(:,c) = mean( Fcell_pix{z}{c}(runBins,:), 2, 'omitnan');
                if somaToggle
                    Froi(runs,z).soma(:,c) = mean( Fsoma_pix{z}{c}(runBins,:), 2, 'omitnan');
                    Froi(runs,z).proc(:,c) = mean( Fproc_pix{z}{c}(runBins,:), 2, 'omitnan');
                end
            end
            % Background subtraction
            Froi(runs,z).cell_sub = Froi(runs,z).cell - Fback{runs,z} + Fback_med(runs,z);
            if somaToggle
                Froi(runs,z).soma_sub = Froi(runs,z).soma - Fback{runs,z} + Fback_med(runs,z);
                Froi(runs,z).proc_sub = Froi(runs,z).proc - Fback{runs,z} + Fback_med(runs,z);
            end

            %{
            figure
            for c = 1:Ncell(z)
                subplot(2,1,1); plot(Froi(runs,z).cell);
                subplot(2,1,2); plot(Froi(runs,z).cell_sub);
            end
            %}

            % Normalization
            % Calculate baseline signals Fo using rolling percentile
            %{
            fprintf('\nCalculating baseline signal');
            Fo = nan(size(F));
            for run = 1:Nruns
                Froi(runs,z).cell_base = MovingPercentile(Froi(runs,z).cell_sub, 10, windowSize, 'pre');
                Fo(runScans{run},:) = MovingPercentile(F(runScans{run},:), basePrct, windowSize, 'pre'); % high-pass filter movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
            end
            Fotot = MovingPercentile(Ftot, basePrct, windowSize, 'pre');
            dFF = (F-Fo)./Fo; % (Ffilt-Fo)./Fo; %dFF(tempCC.PixelIdxList{c},:) = (F-Fo)./Fo;
            dFFtot = (Ftot-Fotot)./Fotot;

            %F = Froi_lp - Fnp_lp + mean(Fnp_lp,1,'omitnan');  % Froi - Fnp + mean(Fnp,1,'omitnan'); % plot( F(:,1) );
            %Ftot = Fall_lp - Fback_lp + mean(Fback_lp,1,'omitnan');
            %}
        end
    end
    % Save the results
    fprintf('\nSaving %s', savePath);
    save(savePath, 'expt', 'projParam', 'refTime', 'Ncell','Froi', 'Trel', 'Fback', 'Fback_med' );
    toc
else
    fprintf('\nLoading %s', savePath);
    load(savePath, 'Froi', 'Trel', 'Fback', 'Ncell');
    toc
end

% Plot the results (optional)
if show
    for z = 1:Nz
        cellColor = distinguishable_colors(Ncell(z));
        fluorescence = figure('WindowState','max');
        for runs = 1:expt.Nruns
            plot( Trel{runs}, Fback{runs,z}, 'k'); hold on;
            for c = 1:Ncell(z)
                plot( Trel{runs}, Froi(runs,z).cell_sub(:,c), 'color', cellColor(c,:)); hold on;
            end
        end
        %{
        h(1) = plot( T, Fback, 'k' ); hold on;
        h(2) = plot( T, Froi.cell, 'b' );
        h(3) = plot( T, Froi.soma, 'r' );
        h(4) = plot( T, Froi.proc, 'g' );
        legend(h, {'Background','Cell','Soma','Processes'});
        %}
        xlim([-Inf,Inf]);
        xlabel('Time (min)'); ylabel('Raw Fluorescence');

        % save the figure
        figPath = sprintf('%s%s_RawFluorescence', expt.dir, expt.name);
        if ~exist(figPath, 'file') || overwrite
            fprintf('\nSaving %s', figPath);
            saveas(fluorescence, figPath)
        end
        pause;
    end
end
end