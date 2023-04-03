function [allBout, allPeri, baseBout, basePeri, stimBout, stimPeri] = FilterPoolPeriBout( locoBout, periBout, mType, boutFilt, varargin ) % allBout, allPeri, 
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
baseBout = cell(1,Nfov); stimBout = cell(1,Nfov); allBout = cell(1,Nfov); 
NgoodBase = zeros(1, Nfov); NgoodStim = zeros(1, Nfov); %totGoodBase = 0; totGoodStim = 0;
for f = flip(1:Nfov) 
    % initialize peri structures
    basePeri(f) = periBout{f}(1); %#ok<*AGROW>
    basePeri(f).Nbout = 0; 
    basePeri(f).velocity = nan(basePeri(f).Nframe, 0);
    basePeri(f).frames = nan(basePeri(f).Nframe, 0);
    basePeri(f).initMax = [];
    basePeri(f).cell = nan(basePeri(f).Nframe, basePeri(f).Ncell, 0);
    basePeri(f).soma = nan(basePeri(f).Nframe, basePeri(f).Ncell, 0);
    basePeri(f).proc = nan(basePeri(f).Nframe, basePeri(f).Ncell, 0);
    basePeri(f).latency = nan(0, basePeri(f).Ncell);
    basePeri(f).peak = nan(0, basePeri(f).Ncell);
    basePeri(f).delta = nan(0, basePeri(f).Ncell);
    basePeri(f).z = nan(0, basePeri(f).Ncell);
    stimPeri(f) = periBout{f}(1);
    stimPeri(f).Nbout = 0; 
    stimPeri(f).velocity = nan(stimPeri(f).Nframe, 0);
    stimPeri(f).frames = nan(stimPeri(f).Nframe, 0);
    stimPeri(f).initMax = [];
    stimPeri(f).cell = nan(stimPeri(f).Nframe, stimPeri(f).Ncell, 0);
    stimPeri(f).soma = nan(stimPeri(f).Nframe, stimPeri(f).Ncell, 0);
    stimPeri(f).proc = nan(stimPeri(f).Nframe, stimPeri(f).Ncell, 0);
    stimPeri(f).latency = nan(0, stimPeri(f).Ncell);
    stimPeri(f).peak = nan(0, stimPeri(f).Ncell);
    stimPeri(f).delta = nan(0, stimPeri(f).Ncell);
    stimPeri(f).z = nan(0, stimPeri(f).Ncell);
    
    % Search each movie for good bouts/peribouts
    % BASELINE MOVIES
    for m = mType(f).base 
        if numel(locoBout{f}{m}) > 0
            % Determine which bouts meet the filter criteria
            durVec = [locoBout{f}{m}.dur]';
            isoMat = vertcat(locoBout{f}{m}.iso); 
            %isoMin = min(isoMat,[],2)';
            bGood = find( isoMat(:,1) >= boutFilt.iso(1) & isoMat(:,2) >= boutFilt.iso(2) & durVec >= boutFilt.dur(1) & durVec <= boutFilt.dur(2) )'; % isoMin >= boutFilt.iso
            NgoodBase(f) = NgoodBase(f) + numel(bGood);
            % Copy the bout and peribout into pooled variables
            baseBout{f} = [baseBout{f}, locoBout{f}{m}(bGood)];
            basePeri(f).frames = [basePeri(f).frames, periBout{f}(m).frames(:,bGood)];
            basePeri(f).velocity = [basePeri(f).velocity, periBout{f}(m).velocity(:,bGood)];
            basePeri(f).initMax = [basePeri(f).initMax, periBout{f}(m).initMax(bGood)];
            basePeri(f).cell = cat(3, basePeri(f).cell, periBout{f}(m).cell(:,:,bGood));
            basePeri(f).soma = cat(3, basePeri(f).soma, periBout{f}(m).soma(:,:,bGood));
            basePeri(f).proc = cat(3, basePeri(f).proc, periBout{f}(m).proc(:,:,bGood));
            basePeri(f).latency = vertcat( basePeri(f).latency, periBout{f}(m).latency.cell(bGood,:) );
            basePeri(f).peak = vertcat( basePeri(f).peak, periBout{f}(m).peak.cell(bGood,:) );
            basePeri(f).delta = vertcat( basePeri(f).delta, periBout{f}(m).delta.cell(bGood,:) );
            basePeri(f).z = vertcat( basePeri(f).z, periBout{f}(m).z.cell(bGood,:) );
            % Plot each bout and peribout (optional)
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
                        h(c) = plot( periBout{f}(m).T, periBout{f}(m).cell(:,c,b) ); hold on;
                        plot( periBout{f}(m).latency.cell(b,c), periBout{f}(m).peak.cell(b,c), 'kx', 'MarkerSize',MS )
                    end
                    h(end+1) = plot( periBout{f}(m).T, mean(periBout{f}(m).cell(:,:,b), 2, 'omitnan'), 'k', 'LineWidth',3 );
                    legend(h, [num2cellStr(1:periBout{f}(m).Ncell), {'Mean'}], 'Location','NorthWest' );
                    xlim([-20,20]);
                    ylabel('dF/F');
                    %dataStr = sprintf('%2.2f  ', periBout{f}(m).z.cell(b,:) ); % periBout{f}(m).delta.cell(b,:)
                    %title( sprintf('Baseline:  [f,m,b] = [%i, %i, %i].  peak z-score = %s', f, m, b, dataStr ) ); % 
                    title( sprintf('Baseline:  [f,m,b] = [%i, %i, %i]. median latency = %2.2f,  median z-score = %2.2f', f, m, b, median(periBout{f}(m).latency.cell(b,:)), median(periBout{f}(m).z.cell(b,:)) ) ); % 
                    pause;
                    clearvars h;
                end
            end
        end
    end
    basePeri(f).Nbout = NgoodBase(f);
    
    % STIMULUS MOVIES
    for m = mType(f).stim 
        if numel(locoBout{f}{m}) > 0
            % Determine which bouts meet the filter criteria
            durVec = [locoBout{f}{m}.dur]';
            isoMat = vertcat(locoBout{f}{m}.iso); 
            %isoMin = min(isoMat,[],2)';
            bGood = find( isoMat(:,1) >= boutFilt.iso(1) & isoMat(:,2) >= boutFilt.iso(2) & durVec >= boutFilt.dur(1) & durVec <= boutFilt.dur(2) )'; % isoMin >= boutFilt.iso
            NgoodStim(f) = NgoodStim(f) + numel(bGood);
            % Copy the bout and peribout into pooled variables
            stimBout{f} = [stimBout{f}, locoBout{f}{m}(bGood)];
            stimPeri(f).frames = [stimPeri(f).frames, periBout{f}(m).frames(:,bGood)];
            stimPeri(f).velocity = [stimPeri(f).velocity, periBout{f}(m).velocity(:,bGood)];
            stimPeri(f).initMax = [stimPeri(f).initMax, periBout{f}(m).initMax(bGood)];
            stimPeri(f).cell = cat(3, stimPeri(f).cell, periBout{f}(m).cell(:,:,bGood));
            stimPeri(f).soma = cat(3, stimPeri(f).soma, periBout{f}(m).soma(:,:,bGood));
            stimPeri(f).proc = cat(3, stimPeri(f).proc, periBout{f}(m).proc(:,:,bGood));
            stimPeri(f).latency = vertcat( stimPeri(f).latency, periBout{f}(m).latency.cell(bGood,:) );
            stimPeri(f).peak = vertcat( stimPeri(f).peak, periBout{f}(m).peak.cell(bGood,:) );
            stimPeri(f).delta = vertcat( stimPeri(f).delta, periBout{f}(m).delta.cell(bGood,:) );
            stimPeri(f).z = vertcat( stimPeri(f).z, periBout{f}(m).z.cell(bGood,:) );
            % Plot each bout and peribout (optional)
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
                        h(c) = plot( periBout{f}(m).T, periBout{f}(m).cell(:,c,b) ); hold on;
                        plot( periBout{f}(m).latency.cell(b,c), periBout{f}(m).peak.cell(b,c), 'kx', 'MarkerSize',MS )
                    end
                    h(end+1) = plot( periBout{f}(m).T, mean(periBout{f}(m).cell(:,:,b), 2, 'omitnan'), 'k', 'LineWidth',3 );
                    legend(h, [num2cellStr(1:periBout{f}(m).Ncell), {'Mean'}], 'Location','NorthWest' );
                    xlim([-20,20]);
                    ylabel('dF/F');
                    %dataStr = sprintf('%2.2f  ', periBout{f}(m).z.cell(b,:) ); % periBout{f}(m).delta.cell(b,:)
                    %title( sprintf('Stimulus:  [f,m,b] = [%i, %i, %i].  peak z-score = %s', f, m, b, dataStr ) ); % 
                    title( sprintf('Stimulus:  [f,m,b] = [%i, %i, %i]. median latency = %2.2f,  median z-score = %2.2f', f, m, b, median(periBout{f}(m).latency.cell(b,:)), median(periBout{f}(m).z.cell(b,:)) ) ); %
                    pause;
                    clearvars h;
                end
            end
        end
    end
    stimPeri(f).Nbout = NgoodStim(f);
    % Merge base/stim pools
    allBout{f} = [baseBout{f}, stimBout{f}]; 
    if NgoodBase(f) > 0 && NgoodStim(f) == 0
        allPeri(f) = basePeri(f);
    elseif NgoodBase(f) == 0 && NgoodStim(f) > 0
        allPeri(f) = stimPeri(f);
    elseif NgoodBase(f) == 0 && NgoodStim(f) == 0
        allPeri(f) = struct('Nbout',0, 'Ncell',NaN, 'Nframe',NaN, 'base', [], 'Nbase', NaN, 'run',[], 'Nrun',NaN, 'Ninit',NaN,...
                'frames',[], 'T',[], 'velocity',[], 'cell',[], 'soma',[], 'proc',[], ...
                'initMax', [], 'latency',[], 'peak',[], 'delta',[], 'z',[] ); 
    else
        allPeri(f) = basePeri(f);
        allPeri(f).Nbout = basePeri(f).Nbout + stimPeri(f).Nbout;
        allPeri(f).velocity = [basePeri(f).velocity, stimPeri(f).velocity];
        allPeri(f).cell = cat(3, basePeri(f).cell, stimPeri(f).cell);
        allPeri(f).soma = cat(3, basePeri(f).soma, stimPeri(f).soma);
        allPeri(f).proc = cat(3, basePeri(f).proc, stimPeri(f).proc);
        allPeri(f).initMax = cat(2, basePeri(f).initMax, stimPeri(f).initMax);
        allPeri(f).latency = cat(1, basePeri(f).latency, stimPeri(f).latency);
        allPeri(f).peak = cat(1, basePeri(f).peak, stimPeri(f).peak);
        allPeri(f).delta = cat(1, basePeri(f).delta, stimPeri(f).delta);
        allPeri(f).z = cat(1, basePeri(f).z, stimPeri(f).z);
    end
end
end