function [ filtPeri, filtBout ] = FilterPoolPeriBout( locoBout, periBout, mType, boutFilt, varargin ) % allBout, allPeri, 
% Search for locomotive bouts with selected properties, and then pool those bouts and peribouts for baseline and stimulus conditions
IP = inputParser;
addRequired( IP, 'locoBout', @iscell )
addRequired( IP, 'periBout', @iscell )
addRequired( IP, 'mType', @isstruct )
addRequired( IP, 'boutFilt', @isstruct )
addParameter( IP, 'show', false, @islogical )
parse( IP, locoBout, periBout, mType, boutFilt, varargin{:} ); 
show = IP.Results.show; 
MS = 8;
if show
    close all; clearvars sp h;
    figure('WindowState','max'); sp(2) = subplot(2,1,2); sp(1) = subplot(2,1,1); 
end
Nfov = numel(mType);
filtBout = cell(1,Nfov); 
typeStr = ''; %#ok<*NASGU>
fprintf('\nFilter criteria: isolation = [%2.1f,  %2.1f], duration = [%2.1f,  %2.1f]', boutFilt.iso(1), boutFilt.iso(2), boutFilt.dur(1), boutFilt.dur(2) );
for f = flip(1:Nfov) 
    goodPeri = find([periBout{f}.Nbout] > 0, 1); % find a non-empty periBout structure
    if isempty(goodPeri), goodPeri = 1; end
    totBout = sum([periBout{f}.Nbout]);
    % initialize filtPeri structures - copy some info from corresponding periBout structure, and initialize empty variables to be filled with the filtered periBouts
    filtPeri(f).Nbout = 0;  %#ok<*AGROW>
    filtPeri(f).Ncell = periBout{f}(goodPeri).Ncell;
    filtPeri(f).Nframe = periBout{f}(goodPeri).Nframe;
    filtPeri(f).frames = periBout{f}(goodPeri).frames;
    filtPeri(f).base = periBout{f}(goodPeri).base;
    filtPeri(f).Nbase = periBout{f}(goodPeri).Nbase;
    filtPeri(f).run = periBout{f}(goodPeri).run;
    filtPeri(f).Nrun = periBout{f}(goodPeri).Nrun;
    filtPeri(f).T = periBout{f}(goodPeri).T;
    filtPeri(f).Ninit = periBout{f}(goodPeri).Ninit;
    filtPeri(f).velocity = nan(filtPeri(f).Nframe, filtPeri(f).Nbout);
    filtPeri(f).frames = nan(filtPeri(f).Nframe, filtPeri(f).Nbout);
    filtPeri(f).dur = nan(filtPeri(f).Nbout,1);
    filtPeri(f).travel = nan(filtPeri(f).Nbout,1); 
    filtPeri(f).maxSpd = nan(filtPeri(f).Nbout,1);
    filtPeri(f).initMax = nan(filtPeri(f).Nbout,1);
    filtPeri(f).type = struct('pre',[], 'Npre',NaN, 'stim',[], 'Nstim',[], 'acute',[], 'Nacute',[]); % specify indices of filtered bouts corresponding to pre-stimulus, post-stimulus, and acute post-stimulus movies
    filtPeri(f).cell.F = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).cell.A = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).cell.dFdt = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).cell.latency = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).cell.peak = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).cell.delta = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).cell.NAD = false(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).cell.trig = false(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).cell.trigFrac = nan(filtPeri(f).Nbout, 1);
    filtPeri(f).soma.F = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).soma.A = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).soma.dFdt = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).soma.latency = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).soma.peak = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).soma.delta = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).soma.NAD = false(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).soma.trig = false(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).soma.trigFrac = nan(filtPeri(f).Nbout, 1);
    filtPeri(f).proc.F = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).proc.A = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).proc.dFdt = nan(filtPeri(f).Nframe, filtPeri(f).Ncell, filtPeri(f).Nbout);
    filtPeri(f).proc.latency = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).proc.peak = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).proc.delta = nan(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).proc.NAD = false(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).proc.trig = false(filtPeri(f).Nbout, filtPeri(f).Ncell);
    filtPeri(f).proc.trigFrac = nan(filtPeri(f).Nbout, 1);
    % Search each movie for good bouts/peribouts
    for m = mType(f).expt %mType(f).base 
        if numel(locoBout{f}{m}) > 0
            % Determine which bouts meet the filter criteria
            durVec = [locoBout{f}{m}.dur]';
            isoMat = vertcat(locoBout{f}{m}.iso); 
            %isoMin = min(isoMat,[],2)';
            bGood = find( isoMat(:,1) >= boutFilt.iso(1) & isoMat(:,2) >= boutFilt.iso(2) & durVec >= boutFilt.dur(1) & durVec <= boutFilt.dur(2) )'; % isoMin >= boutFilt.iso
            Ngood = numel(bGood);
            if Ngood > 0
                % Enter the bouts' indices into the appropriate movie-type substructure
                bChunk = (1:Ngood) + filtPeri(f).Nbout;
                if ismember(m, mType(f).base)
                    filtPeri(f).type.pre = [filtPeri(f).type.pre, bChunk];
                    typeStr = 'Pre-stimulus';
                elseif ismember(m, mType(f).stim)
                    filtPeri(f).type.stim = [filtPeri(f).type.stim, bChunk];
                    typeStr = 'Post-stimulus';
                else
                    error('Undefined movie type');
                end
                if numel(mType(f).acute) == 2 && m == mType(f).acute(2) 
                    filtPeri(f).type.acute = [filtPeri(f).type.acute, bChunk];
                    typeStr = 'Acute post-stimulus';
                end
                % Copy the bout and peribout into pooled variables
                filtBout{f} = [filtBout{f}, locoBout{f}{m}(bGood)];
                filtPeri(f).Nbout = filtPeri(f).Nbout + Ngood;
                filtPeri(f).frames = [filtPeri(f).frames, periBout{f}(m).frames(:,bGood)];
                filtPeri(f).velocity = [filtPeri(f).velocity, periBout{f}(m).velocity(:,bGood)];
                filtPeri(f).dur = [filtPeri(f).dur; periBout{f}(m).dur(bGood)];
                filtPeri(f).travel = [filtPeri(f).travel; periBout{f}(m).travel(bGood)];
                filtPeri(f).maxSpd = [filtPeri(f).maxSpd; periBout{f}(m).maxSpd(bGood)];
                filtPeri(f).initMax = [filtPeri(f).initMax; periBout{f}(m).initMax(bGood)];
                filtPeri(f).cell.F = cat(3, filtPeri(f).cell.F, periBout{f}(m).cell.F(:,:,bGood));
                filtPeri(f).cell.A = cat(3, filtPeri(f).cell.A, periBout{f}(m).cell.A(:,:,bGood));
                filtPeri(f).cell.latency = vertcat( filtPeri(f).cell.latency, periBout{f}(m).cell.latency(bGood,:) );
                filtPeri(f).cell.peak = vertcat( filtPeri(f).cell.peak, periBout{f}(m).cell.peak(bGood,:) );
                filtPeri(f).cell.delta = vertcat( filtPeri(f).cell.delta, periBout{f}(m).cell.delta(bGood,:) );
                filtPeri(f).cell.NAD = vertcat( filtPeri(f).cell.NAD, periBout{f}(m).cell.NAD(bGood,:) );
                filtPeri(f).cell.trig = vertcat( filtPeri(f).cell.trig, periBout{f}(m).cell.trig(bGood,:) );
                filtPeri(f).cell.trigFrac = vertcat( filtPeri(f).cell.trigFrac, periBout{f}(m).cell.trigFrac(bGood) );
                filtPeri(f).soma.F = cat(3, filtPeri(f).soma.F, periBout{f}(m).soma.F(:,:,bGood));
                filtPeri(f).soma.A = cat(3, filtPeri(f).soma.A, periBout{f}(m).soma.A(:,:,bGood));
                filtPeri(f).soma.latency = vertcat( filtPeri(f).soma.latency, periBout{f}(m).soma.latency(bGood,:) );
                filtPeri(f).soma.peak = vertcat( filtPeri(f).soma.peak, periBout{f}(m).soma.peak(bGood,:) );
                filtPeri(f).soma.delta = vertcat( filtPeri(f).soma.delta, periBout{f}(m).soma.delta(bGood,:) );
                filtPeri(f).soma.NAD = vertcat( filtPeri(f).soma.NAD, periBout{f}(m).soma.NAD(bGood,:) );
                filtPeri(f).soma.trig = vertcat( filtPeri(f).soma.trig, periBout{f}(m).soma.trig(bGood,:) );
                filtPeri(f).soma.trigFrac = vertcat( filtPeri(f).soma.trigFrac, periBout{f}(m).soma.trigFrac(bGood) );
                filtPeri(f).proc.F = cat(3, filtPeri(f).proc.F, periBout{f}(m).proc.F(:,:,bGood));
                filtPeri(f).proc.A = cat(3, filtPeri(f).proc.A, periBout{f}(m).proc.A(:,:,bGood));
                filtPeri(f).proc.latency = vertcat( filtPeri(f).proc.latency, periBout{f}(m).proc.latency(bGood,:) );
                filtPeri(f).proc.peak = vertcat( filtPeri(f).proc.peak, periBout{f}(m).proc.peak(bGood,:) );
                filtPeri(f).proc.delta = vertcat( filtPeri(f).proc.delta, periBout{f}(m).proc.delta(bGood,:) );
                filtPeri(f).proc.NAD = vertcat( filtPeri(f).proc.NAD, periBout{f}(m).proc.NAD(bGood,:) );
                filtPeri(f).proc.trig = vertcat( filtPeri(f).proc.trig, periBout{f}(m).proc.trig(bGood,:) );
                filtPeri(f).proc.trigFrac = vertcat( filtPeri(f).proc.trigFrac, periBout{f}(m).proc.trigFrac(bGood) );
                if show && ~isempty( bGood )
                    for b = bGood
                        subplot(sp(2)); cla;
                        plot( periBout{f}(m).T, periBout{f}(m).velocity(:,b) ); hold on; % locoBout{f}{m}.T locoBout{f}{m}.velocity(:,b)
                        xlim([-20,20]);
                        ylabel('Velocity (cm/s)'); xlabel('Peri-bout Time (s)');
                        title( sprintf('iso = [%2.1f, %2.1f], duration = %2.2f s', locoBout{f}{m}(b).iso(1), locoBout{f}{m}(b).iso(2), locoBout{f}{m}(b).dur ) )
                        subplot(sp(1)); cla;
                        %plot( periBout{f}(m).T, periBout{f}(m).cell(:,:,b) ); hold on;
                        for c = flip(1:periBout{f}(m).Ncell)
                            h(c) = plot( periBout{f}(m).T, periBout{f}(m).cell.A(:,c,b) ); hold on;
                            plot( periBout{f}(m).cell.latency(b,c), periBout{f}(m).cell.peak(b,c), 'kx', 'MarkerSize',MS )
                        end
                        h(end+1) = plot( periBout{f}(m).T, mean(periBout{f}(m).cell.A(:,:,b), 2, 'omitnan'), 'k', 'LineWidth',3 );
                        legend(h, [num2cellStr(1:periBout{f}(m).Ncell), {'Mean'}], 'Location','NorthWest' );
                        xlim([-20,20]);
                        ylabel('Cell Activity');
                        %dataStr = sprintf('%2.2f  ', periBout{f}(m).z.cell(b,:) ); % periBout{f}(m).cell.delta(b,:)
                        title( sprintf('%s [f,m,b] = [%i, %i, %i]. median latency = %2.2f,  median z-score = %2.2f', typeStr, f, m, b, median(periBout{f}(m).cell.latency(b,:)), median(periBout{f}(m).z.cell(b,:)) ) ); % 
                        pause;
                        clearvars h;
                    end
                end
            end
        end
    end
    filtPeri(f).type.Npre = numel( filtPeri(f).type.pre ); filtPeri(f).type.Nstim = numel( filtPeri(f).type.stim ); filtPeri(f).type.Nacute = numel( filtPeri(f).type.acute ); 
    dT = median( diff(filtPeri(f).T) );
    filtPeri(f).cell.dFdt = diff( filtPeri(f).cell.F, 1, 1 )/dT; filtPeri(f).soma.dFdt = diff( filtPeri(f).soma.F, 1, 1 )/dT; filtPeri(f).proc.dFdt = diff( filtPeri(f).proc.F, 1, 1 )/dT;
    fprintf('\n  f = %i:  %i  of  %i  periBouts met filtering criteria', f, filtPeri(f).Nbout, totBout );
end
fprintf('\n\n');
end