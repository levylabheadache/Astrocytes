clear; clc; close all;
% Load the master spreadsheet
dataDir = 'D:\2photon\'; %'D:\Andy\AstroCa\';
dataTablePath = 'D:\MATLAB\Astrocytes\AstrocyteData.xlsx';
dataTable = readtable(dataTablePath, 'sheet', 'Experiments');
GetRunNumber = @(x)(str2double(x(strfind(x, 'run')+3:end)));
% Find the specific experiments needed
%setGroup = {'Awake/Baseline', 'Awake/Saline', 'Awake/CNO/Ctrl', 'Awake/CLZ', 'Low CNO', 'Medium CNO', 'High CNO', 'Awake/C21', }; %{'Medium CNO', 'Medium CNO IP'}; % {'Iso/NTG', 'Awake/NTG'}; %  'Low CNO', 'Medium CNO', 'High CNO', 'Awake/NTG'}; %   {'C21'}; %  'Low CNO', 'Medium CNO', 'High CNO'    'Iso/Saline','Iso/Medium CNO','Iso/C21', 'Iso/Saline','Iso/Medium CNO', 'Saline', 'Medium CNO'
%Ngroup = numel(setGroup);
%for g = flip(1:Ngroup), groupRows{g} = find( strcmp( exptTable.Group, setGroup{g} ) )'; end
%NgroupExpt = cellfun(@numel, groupRows);
%allRows = [groupRows{:}];
Nexpt = size(dataTable, 1); %numel( allRows );  % sum( NgroupExpt ); % size(exptTable,1);
% registration parameters
regParams.refChan = 'red'; %green';  % 
regParams.refRun = 3;
regParams.refScan = 300:700; %1508:1655; %32280:43420; % 36910:43430; % 100,100,20,20.
regParams.chunkSize = 1000; %
regParams.histmatch = true; % false; %
regParams.avgT = 15;
regParams.avgTsigma = 5;
regParams.binXY = 2;
regParams.binT = 1;
regParams.prereg = false;
regParams.highpass = 0;
regParams.lowpass = 0;
regParams.medFilter = [0,0,0];
regParams.minInt = 1500;
regParams.edges = [100,100,40,40]; %[120,100,30,30]; % [80,90,20,20]; %[60,30,20,20]; % [80,90,20,20]; % [20,60,20,20];
regParams.name = '';
regParams.method = 'translation'; % 'affine';

% Gather the relevant data from each experiment
for x = 45%1:Nexpt  %1:Nexpt
    % Parse data table
    expt(x).mouse = dataTable.Mouse{x}; %dataTable{x,dataCol.mouse};
    expt(x).date = dataTable.Date(x); %dataTable{x,dataCol.date};  
    if isnumeric(expt(x).date), expt(x).date = num2str(expt(x).date); end
    expt(x).fov = str2double(dataTable.FOV(x)) ; %dataTable{x,dataCol.FOV};
    expt(x).dir = sprintf('%s%s\\%s_FOV%i\\', dataDir, expt(x).mouse, expt(x).date, expt(x).fov); % sprintf('%s%s\\%s_FOV%i_%s\\', dataDir, expt(x).mouse, expt(x).date, expt(x).fov, expt(x).mouse);
    expt(x).name = sprintf('%s_%s_FOV%i', expt(x).mouse, expt(x).date, expt(x).fov);
    runFolders = FileFinder(expt(x).dir, 'contains','run', 'type',0);
    [expt(x).runs, sortInd] = sort(cellfun(GetRunNumber, runFolders, 'UniformOutput',true)', 'ascend');
    expt(x).Nruns = numel(expt(x).runs); %expt(x).Nroi = NaN; expt(x).Naxon = NaN;
    expt(x).drug = dataTable.Drug{x};   
    expt(x).dose = dataTable.Dose_mg_kg_(x); 

    % Get and concatenate run-level metadata and times
    runInfo{x} = GetRunInfo(expt(x), dataDir);
    for run = 1:expt(x).Nruns, expt(x).runMovies(run) = str2double(runInfo{x}(run).fileName(end-2:end));  end
    expt(x).baseRuns = find(expt(x).runMovies <= dataTable.LastBaseMovie(x));
    expt(x).stimRuns = find(expt(x).runMovies >= dataTable.FirstStimMovie(x));
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x});
    expt(x).Nrow = runInfo{x}(1).sz(1); expt(x).Ncol = runInfo{x}(1).sz(2); expt(x).Nplane = runInfo{x}(1).otlevels; expt(x).Nchan = runInfo{x}(1).nchan;
    expt(x).Nscan = floor([runInfo{x}.nframes]/expt(x).Nplane); expt(x).totScan = sum(expt(x).Nscan); expt(x).totFrame = sum([runInfo{x}.totFrame]);
    expt(x).scanLims = [0, cumsum(expt(x).Nscan)];
    expt(x).frameRate = runInfo{x}(1).framerate; expt(x).scanRate = runInfo{x}(1).framerate/expt(x).Nplane;
    catInfo{x} = ConcatenateRunInfo(expt(x), runInfo{x}, 'suffix','sbxcat', 'overwrite',false); % Get concatenated metadata
    expt(x).zoom = str2double(catInfo{x}.config.magnification_list(catInfo{x}.config.magnification,:)); % fprintf('\nMagnification = %2.1f X', expt(x).zoom);
    expt(x).umPerPixel = (1/0.53)/expt(x).zoom;
    expt(x).sbx = strcat(expt(x).dir, expt(x).name, '.sbxcat'); % sbx_interp

    % Perform registration on individual runs, and concatenate into a single sbx file
    if ~exist(expt(x).sbx, 'file')
        for t = flip(expt(x).runs)
            RegisterCat3D( runInfo{x}(t), 'overwrite',false, 'writeZ',false, 'Zint',3, 'chunk',20, 'edge',[90,140,60,20], 'refChan','green', 'preaff',true, 'fix',true);  % 1500
        end
        close all;
        ConcatenateExptRuns(expt(x), runInfo{x}, catInfo{x}, 'refRun',1, 'setEdge',[100,100,40,40]);
    end

    % Run registration on concatenated data
    RegisterCat3D( catInfo{x}, regParams, 'overwrite',false, 'writeZ',false, 'chunk',30, 'preaff',false); % expt(x).mouse, expt(x).date, expt(x).runs , 'cat',false
    [~,regPath] = FileFinder(expt(x).dir, 'type','sbxreg');
    if ~isempty(regPath), expt(x).sbx = regPath{1}; end %expt(x).sbx = strcat(expt(x).dir, expt(x).name, '.sbxreg');
     
    % Generate downsampled movies for each run from the concatenated data
    projParam.dir = strcat(expt(x).dir, 'Projections\');
    projParam.bin = 8;
    projParam.scale = 4;
    %GetEdge()
    projParam.edge = [60,60,40,20]; % [left, right, top, bottom]
    mkdir( projParam.dir)
    projParamPath = sprintf('%s%s_projParam.mat',projParam.dir, expt(x).name );
    greenCatProjPath = sprintf('%s%s_greenProj_cat.tif', projParam.dir, expt(x).name );
    WriteSbxPlaneTif(expt(x).sbx, catInfo{x}, 1, 'chan','green', 'firstScan',8, ...
        'edge',projParam.edge, 'scale',projParam.scale, 'binT',projParam.bin, 'RGB',false, 'dir',[projParam.dir,'Green\'], 'name',sprintf('%s_cat',expt(x).name), 'overwrite',false );
    if ~exist(greenCatProjPath, 'file')
        save(projParamPath, 'projParam');
        greenRunProj = cell(1,expt(x).Nruns); %redRunProj = cell(1,expt(x).Nruns);
        for run = 1:expt(x).Nruns
            if expt(x).Nplane == 1
                greenRunProj{run} = WriteSbxPlaneTif(expt(x).sbx, catInfo{x}, 1, 'chan','green', 'firstScan',expt(x).scanLims(run)+projParam.bin+1, 'Nscan',projParam.bin*floor((expt(x).Nscan(run)-projParam.bin)/projParam.bin), ...
                    'edge',projParam.edge, 'scale',projParam.scale, 'binT',projParam.bin, 'RGB',false, 'dir',[projParam.dir,'Green\'], 'name',sprintf('%s_run%i',expt(x).name, run), 'overwrite',false );
                %[~,runProj{run}] = WriteSbxPlaneTif(expt(x).sbx, catInfo{x}, 1, 'chan','red', 'firstScan',expt(x).scanLims(run)+projParam.bin+1, 'Nscan',projParam.bin*floor((expt(x).Nscan(run)-projParam.bin)/projParam.bin), ...
                %    'edge',projParam.edge, 'scale',projParam.scale, 'binT',projParam.bin, 'RGB',false, 'dir',[projParam.dir,'Red\'], 'name',sprintf('%s_run%i',expt(x).name, run), 'overwrite',false );
            end
        end
        % Concatenate projections
        %runProj = vertcat(greenRunProj{expt(x).runs});
        fprintf('Writing %s', greenCatProjPath)
        WriteTiff(cat(3, greenRunProj{:}), greenCatProjPath );
        %WriteTiff(cat(3, redRunProj{:}), sprintf('%s%s_redProj_cat.tif', projParam.dir, expt(x).name ) );
    else
        load(projParamPath)
    end
    
    % Draw astrocyte ROIs
    %[backROI{x}, cellROI{x}, somaROI{x}, procROI{x}, Ncell{x}, proj{x}] = MakeAstroROI( expt(x), 'overwrite',false );
    %[cellEvt{x}, aquaOpts{x}, aquaResults{x}] = AquaBatchRuns(expt(x), projParam.bin/expt(x).scanRate, expt(x).umPerPixel*projParam.scale, 0);
end

%%
timeStamp = datetime;
NgroupFOV = nan(1,Ngroup); NgroupAstro = nan(1,Ngroup);
uniqueMice = cell(1,Ngroup); NuniqueMice = nan(1,Ngroup); NgroupM = nan(1,Ngroup); NgroupF = nan(1,Ngroup);
for g = 1:Ngroup
    NgroupFOV(g) = sum(Nfov(xGroup{g}));
    NgroupAstro(g) =  sum(Nastro(xGroup{g}));
    [uniqueMice{g}, uniqueInd ] = unique( mouse(xGroup{g}) );
    NuniqueMice(g)  = numel( uniqueInd);
    uniqueSex = sex(xGroup{g}(uniqueInd));
    NgroupM(g) = sum(strcmpi( uniqueSex, 'M' ));
    NgroupF(g) = sum(strcmpi( uniqueSex, 'F' ));
    fprintf('\n%s:  %i experiments, %i FoV, %i mice (%i M / %i F), %i astrocytes' , setGroup{g}, NgroupExpt(g), NgroupFOV(g), NuniqueMice(g), NgroupM(g), NgroupF(g), NgroupAstro(g) );
end

censusTable = table( setGroup', NgroupExpt', NgroupFOV', NgroupAstro', NuniqueMice', NgroupM', NgroupF', repmat(timeStamp,Ngroup,1), 'VariableNames',{'Group','Nexpt', 'Nfov', 'Nastro','Nmice','Nmale','Nfemale','TimeStamp'},'RowNames',setGroup ); 
writetable(censusTable, 'D:\2photon\AstroDataCensus.xls');

%% Calculate summary statistics from the relevant movies
Fmed = cell(1,Nexpt); FbackMed = cell(1,Nexpt); Fnorm = cell(1,Nexpt); 
lastObs = cell(1,Nexpt);
for x = 1:Nexpt
    for f = 1:Nfov(x)
        % Calculate median values
        FbackMed{x}{f} = cellfun(@median, Fback{x}{f} )';
        tempMed = cellfun( @median, {Fastro{x}{f}.cell}, 'UniformOutput',false );
        Fmed{x}{f} = vertcat(tempMed{:});
        Fsub = Fmed{x}{f} - repmat( FbackMed{x}{f}, 1, Ncell{x}(f) );
        FsubBase = mean( Fsub(mBase{x}{f},:), 1 );
        Fnorm{x}{f} = Fsub./repmat( FsubBase, size(Fsub,1), 1 );
         
        lastObs{x}(f) = Tmed{x}{f}(end); % mStim{x}{f}(end)
    end
end

%% time course: full fluor data plus locomotion
Tshow = cell(1,Nexpt); sp = cell(1,Nexpt);
close all; 
opt = {[0.06,0.03], [0.08,0.04], [0.07,0.07]};  % {[vert, horz], [bottom, top], [left, right] } 
for x = 1:Nexpt
    % Set up the figure
    figure('units','normalized','WindowState','max', 'color','w');
    for f = 1:Nfov(x)+1, sp{x}(f) = subtightplot(Nfov(x)+1, 1, f, opt{:}); end
    % Plot the data
    for f = flip(1:Nfov(x))
        % Fluor
        Tshow{x}{f} = vertcat(Tmin{x}{f}{mExpt{x}{f}});
        subplot(sp{x}(f)) %subtightplot(Nfov(x)+1, 1, f, opt{:});
        plot( Tshow{x}{f}, vertcat(Fback{x}{f}{mExpt{x}{f}}), 'color',0.5*[1,1,1] ); hold on;
        plot( Tshow{x}{f}, vertcat(Fastro{x}{f}(mExpt{x}{f}).cell) );
        ylabel('Raw Fluorescence'); title( sprintf('FOV %i', fov{x}(f) ) );
        % Locomotion
        subplot(sp{x}(Nfov(x)+1))
        if ~any(cellfun( @isempty, velocity{x}{f} ))
            for m = mExpt{x}{f}, plot( Tmin{x}{f}{m}, velocity{x}{f}{m}, 'k' ); hold on;  end
        end
        %pause;
    end
    xlabel( 'Time After Injection (min)' ); ylabel( 'Velocity (cm/s)' ); 
    
    linkaxes(sp{x}, 'x'); axis tight;
    subplot(sp{x}(1));  
    xlim( [min(cellfun(@min, Tshow{x} )), max(cellfun(@max, Tshow{x} ))]);
    title( sprintf('%s %s', mouse{x}, exptDate{x} ) ); % sprintf('%s %s: %2.2f mg/kg %s.  FOV %i', mouse{x}, exptDate{x}, dose{x}, drug{x}, fov{x}(f) )
    pause;
end


%% time course: full fluor data only
close all; clearvars sp h;
for x = 1:Nexpt
    figure('units','normalized','WindowState','max', 'color','w');
    fovColor = cell(1,Nfov(x));
    fullColor = distinguishable_colors( sum(Ncell{x}) );
    for f = 1:Nfov(x)
        fovColor{f} = fullColor(1:Ncell{x}(f),:);
        fullColor(1:Ncell{x}(f),:) = [];
        for m = [mBase{x}{f}, mStim{x}{f}]
            plot( Tmin{x}{f}{m}, Fback{x}{f}{m}, 'k' ); hold on;
            for c = 1:Ncell{x}(f)
                plot( Tmin{x}{f}{m}, Fastro{x}{f}(m).cell(:,c), 'color',fovColor{f}(c,:) ); hold on;
            end
            
        end
        %pause;
    end
    xlabel( sprintf('Time After %s Injection (min)', drug{x}) ); ylabel('Raw Fluorescence');
    title( sprintf('%s %s', mouse{x}, exptDate{x} ) );
    axis tight;
    pause;
end

%% time course: normalized median fluorescence
close all; clearvars sp h;
for x = Nexpt
    figure('units','normalized','WindowState','max', 'color','w');
    fovColor = cell(1,Nfov(x));
    fullColor = distinguishable_colors( sum(Ncell{x}) );
    for f = 1:Nfov(x)
        fovColor{f} = fullColor(1:Ncell{x}(f),:);
        fullColor(1:Ncell{x}(f),:) = [];
        for c = 1:Ncell{x}(f)
            plot(Tmed{x}{f}(mExpt{x}{f}), Fnorm{x}{f}(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:)); hold on;
            plot(Tmed{x}{f}(mExpt{x}{f}), Fnorm{x}{f}(mExpt{x}{f},c), 'color',fovColor{f}(c,:));
        end
    end
    xlabel('Time Post Injection (min)'); ylabel('Median Fluorescence');
    %title( sprintf('%s %s  (%s %1.1f mg/kg)', mouse{x}, exptDate{x}, drug{x}, exptTable.Dose_mg_kg_(x) ) );
    axis tight;
    pause;
end

%% time course: normalized median fluorescence by drug
exptColor = distinguishable_colors(Nexpt);
alphaVal = 0.4;
close all; clearvars sp h;
opt = {[0.07,0.05], [0.1,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
figure('units','normalized','WindowState','max', 'color','w');
for g = 1:numel(drugTypes)
    sp(g) = subtightplot(NdrugTypes,1,g,opt{:});
    for x = xDrug{g} %10% 1:Nexpt
        for f = 1:Nfov(x)
            for c = 1:Ncell{x}(f)
                plot(Tmed{x}{f}(mExpt{x}{f}), Fnorm{x}{f}(mExpt{x}{f},c), '.', 'color',[exptColor(x,:),1]); hold on;
                plot(Tmed{x}{f}(mExpt{x}{f}), Fnorm{x}{f}(mExpt{x}{f},c), 'color',[exptColor(x,:),alphaVal]);
            end
        end
    end
    ylabel('Median Normalized Fluorescence');
    title( sprintf('%s', drugTypes{g}) );
end
xlabel('Time Post-Injection (min)'); 
linkaxes(sp,'xy');

%% time course: soma/cell fluorescence
Tmed = cell(1,Nexpt); somaCellRatio = cell(1,Nexpt);
close all; clearvars sp h;
for x = 10 %1:Nexpt
    figure('units','normalized','WindowState','max', 'color','w');
    fovColor = cell(1,Nfov(x));
    fullColor = distinguishable_colors( sum(Ncell{x}) );
    for f = 1:Nfov(x)
        fovColor{f} = fullColor(1:Ncell{x}(f),:);
        fullColor(1:Ncell{x}(f),:) = [];
        mExpt{x}{f} = [mBase{x}{f}, mStim{x}{f}];
        % Calculate median values
        Tmed{x}{f} = cellfun( @median, Tmin{x}{f} )';
        tempCellMed = cellfun( @median, {Fastro{x}{f}.cell}, 'UniformOutput',false );
        tempSomaMed = cellfun( @median, {Fastro{x}{f}.soma}, 'UniformOutput',false );
        tempRatio = cellfun( @(X,Y)X./Y, tempSomaMed, tempCellMed, 'UniformOutput',false );
        somaCellRatio{x}{f} = vertcat( tempRatio{:} );
        for c = 1:Ncell{x}(f)
            plot(Tmed{x}{f}(mExpt{x}{f}), somaCellRatio{x}{f}(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:)); hold on;
            plot(Tmed{x}{f}(mExpt{x}{f}), somaCellRatio{x}{f}(mExpt{x}{f},c), 'color',fovColor{f}(c,:)); hold on;
        end
        %pause;
    end
    xlabel('Time Post-Injection (min)'); ylabel('Soma/Cell Ratio');
    title( sprintf('CNO Data: %s %s  (%s)', mouse{x}, exptDate{x}, exptTable.StimulusDose{x} ) );
    pause;
end

%% time course: low-freq power
Tmed = cell(1,Nexpt); somaCellRatio = cell(1,Nexpt);
close all; clearvars sp h;
for x = 1:Nexpt
    figure('units','normalized','WindowState','max', 'color','w');
    fovColor = cell(1,Nfov(x));
    fullColor = distinguishable_colors( sum(Ncell{x}) );
    for f = 1:Nfov(x)
        fovColor{f} = fullColor(1:Ncell{x}(f),:);
        fullColor(1:Ncell{x}(f),:) = [];
        mExpt{x}{f} = [mBase{x}{f}, mStim{x}{f}];
        % Calculate median values
        Tmed{x}{f} = cellfun( @median, Tmin{x}{f} )';
        for c = 1:Ncell{x}(f)
            plot(Tmed{x}{f}(mExpt{x}{f}), LFfrac{x}{f}(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:)); hold on;
            plot(Tmed{x}{f}(mExpt{x}{f}), LFfrac{x}{f}(mExpt{x}{f},c), 'color',fovColor{f}(c,:)); hold on;
        end
        %pause;
    end
    xlabel('Time Post-Injection (min)'); ylabel('Low Frequency Fraction');
    title( sprintf('CNO Data: %s %s  (%s)', mouse{x}, exptDate{x}, exptTable.StimulusDose{x} ) );
    pause;
end

%% before/after: low-freq power

close all; clearvars sp h;
%opt = {[0.12,0.05], [0.1,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
figure('WindowState','normal', 'color','w');
shiftX = 0.4;
for x = 1:Nexpt
    for f = 1:Nfov(x)
        cla;
        cellColor = distinguishable_colors(Ncell{x}(f));
        % Gather baseline datapoints
        LFbase = LFfrac{x}{f}( mBase{x}{f}, : );
        LFstim = LFfrac{x}{f}( mStim{x}{f}, : );
        for c = 1:Ncell{x}(f)
            plot([c,c+shiftX],[mean(LFbase(:,c),1), mean(LFstim(:,c),1)], 'Color',cellColor(c,:) ); hold on;
            errorbar( c, mean(LFbase(:,c),1), SEM(LFbase(:,c)), 'Color',cellColor(c,:), 'CapSize',0, 'LineWidth',0.1 );
            errorbar( c+shiftX, mean(LFstim(:,c),1), SEM(LFstim(:,c)), 'Color',cellColor(c,:), 'CapSize',0, 'LineWidth',0.1 );
            plot( c, LFbase(:,c), 'o', 'Color',cellColor(c,:) ); hold on;
            plot( c+shiftX, LFstim(:,c), 'x', 'Color',cellColor(c,:) );
        end
        set(gca, 'Xtick',1:Ncell{x}(f), 'TickDir','out','box','off' );
        xlim([0.5, Ncell{x}(f)+0.5]);
        xlabel('Cell #'); ylabel('Low Frequency Power Fraction'); 
        title( sprintf('CNO Data: %s %s FOV %i (%s)', mouse{x}, exptDate{x}, fov{x}(f), dose{x} ) );
        pause; 
    end
end


%% spectrograms for each cell/movie
close all; clearvars sp h;
opt = {[0.04,0.05], [0.1,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
figure('units','normalized','WindowState','max', 'color','w');
for x = 10 %1:Nexpt
    for f = flip(1:Nfov(x))
        Nshow = numel(mBase{x}{f}) + numel(mStim{x}{f}); 
        for c = 1:Ncell{x}(f)
            m = 0;
            for m = [mBase{x}{f}, mStim{x}{f}]
                m = m+1;
                subtightplot(Nshow,1,m, opt{:}); cla;
                plot_matrix( Pgram{x}{f}{m}(:,:,c), Tgram{x}{f}{m}, fGram{x}{f}{m} );  colorbar('off');
                ylabel('Frequency (Hz)'); title( sprintf('[x,f,c,m] = [%i, %i, %i, %i]', x, f, c, m) );
            end
            xlabel('Time (s)');
            pause;
        end
    end
    
end

%% Time course of Ca events (AQuA)
close all; clearvars sp h;
for x = 1:Nexpt
    figure('units','normalized','WindowState','max', 'color','w');
    fovColor = cell(1,Nfov(x));
    fullColor = distinguishable_colors( sum(Ncell{x}) );
    for f = 1:Nfov(x)
        fovColor{f} = fullColor(1:Ncell{x}(f),:);
        fullColor(1:Ncell{x}(f),:) = [];
        sp(1) = subplot(4,1,1);
        for c = 1:Ncell{x}(f) 
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.rate(mExpt{x}{f},c), 'color',fovColor{f}(c,:) ); hold on;
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.rate(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:) );
        end
        ylabel('Event Rate (per minute)'); title( sprintf('%s %s', mouse{x}, exptDate{x} ) );
        
        sp(2) = subplot(4,1,2);
        for c = 1:Ncell{x}(f) 
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.area(mExpt{x}{f},c), 'color',fovColor{f}(c,:) ); hold on;
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.area(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:) );
        end
        ylabel('Mean Event Area (um^2)');
        
        sp(3) = subplot(4,1,3);
        for c = 1:Ncell{x}(f) 
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.dur(mExpt{x}{f},c), 'color',fovColor{f}(c,:) ); hold on;
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.dur(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:) );
        end
        ylabel('Mean Duration (s)'); 
        
        sp(4) = subplot(4,1,4);
        for c = 1:Ncell{x}(f) 
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.dFF(mExpt{x}{f},c), 'color',fovColor{f}(c,:) ); hold on;
            plot( Tmed{x}{f}(mExpt{x}{f}), cellEvtSumm{x}{f}.dFF(mExpt{x}{f},c), '.', 'color',fovColor{f}(c,:) );
        end
        ylabel('Mean dF/Fo'); 
        xlabel( sprintf('Time After Injection (min)') ); 
    end
    linkaxes(sp,'x'); axis tight;
end

%% Selected example kymographs pre-injection/30 mins/60 mins
saveDir = 'D:\MATLAB\Figures\';
RGBOpt = struct('overwrite',true, 'message',true, 'append',false, 'big',false, 'color',true );
opt = {[0.06,0.01], [0.1,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
colorLims = [0.5,8];
frameLim = [1, round(4*60*15.49/3)];

close all; % clearvars sp h;
ExampleKymographs = figure('units','normalized','WindowState','max', 'color','w');
% ISO/SALINE
X = 2; F = 2; M = [3,5,8]; % MTG003 200122 FOV3
baseFluor = prctile( Fsub{X}{F}(M(1)).cell, 10 );
% Pre-injection
subtightplot(4,3,1,opt{:});
imagesc( (Fsub{X}{F}(M(1)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'ISO_SALINE_PRE.tif'], '-dtiff', '-r200' ); 
title('PRE-INJECTION'); ylabel('ISO/SALINE', 'FontWeight','bold');

% ~30 mins post
subtightplot(4,3,2,opt{:});
imagesc( (Fsub{X}{F}(M(2)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'ISO_SALINE_30MIN.tif'], '-dtiff', '-r200' ); 
title('30 MINS POST-INJECTION');

% ~60 mins post
subtightplot(4,3,3,opt{:});
imagesc( (Fsub{X}{F}(M(3)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
colorbar;
%print( ExampleKymographs, [saveDir,'ISO_SALINE_60MIN.tif'], '-dtiff', '-r200' ); 
title('60 MINS POST-INJECTION');
%cb.Label.String = 'Baseline-Normalized Mean Fluorescence';
%cb.FontSize = FS;

% ISO/CNO
X = 4; F = 1; M = [3,6,8]; % AB26 200108 FOV1
baseFluor = prctile( Fsub{X}{F}(M(1)).cell, 10 );
% Pre-injection
subtightplot(4,3,4,opt{:});
imagesc( (Fsub{X}{F}(M(1)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'ISO_CNO_PRE.tif'], '-dtiff', '-r200' ); 
ylabel('ISO/CNO', 'FontWeight','bold');
% ~30 mins post
subtightplot(4,3,5,opt{:});
imagesc( (Fsub{X}{F}(M(2)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'ISO_CNO_30MIN.tif'], '-dtiff', '-r200' ); 
% ~60 mins post
subtightplot(4,3,6,opt{:});
imagesc( (Fsub{X}{F}(M(3)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
colorbar;
%print( ExampleKymographs, [saveDir,'ISO_CNO_60MIN.tif'], '-dtiff', '-r200' ); 


% AWAKE/SALINE
X = 8; F = 1; M = [3,6,8]; % AB23 191122 FOV1
baseFluor = prctile( Fsub{X}{F}(M(1)).cell, 10 );
% Pre-injection
subtightplot(4,3,7,opt{:});
imagesc( (Fsub{X}{F}(M(1)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'AWAKE_SALINE_PRE.tif'], '-dtiff', '-r200' );
ylabel('AWAKE/SALINE', 'FontWeight','bold');

% ~30 mins post
subtightplot(4,3,8,opt{:});
imagesc( (Fsub{X}{F}(M(2)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'AWAKE_SALINE_30MIN.tif'], '-dtiff', '-r200' );

% ~60 mins post
subtightplot(4,3,9,opt{:});
imagesc( (Fsub{X}{F}(M(3)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
colorbar;
%print( ExampleKymographs, [saveDir,'AWAKE_SALINE_60MIN.tif'], '-dtiff', '-r200' );

% AWAKE/CNO
%X = 13; F = 1; M = [2,5,6];
X = 13; F = 2; M = [2,4,5]; % AB26 191223 FOV2
baseFluor = prctile( Fsub{X}{F}(M(1)).cell, 10 );

% Pre-injection
subtightplot(4,3,10,opt{:});
imagesc( (Fsub{X}{F}(M(1)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'AWAKE_CNO_PRE.tif'], '-dtiff', '-r200' );
ylabel('AWAKE/CNO', 'FontWeight','bold');

% ~30 mins post
subtightplot(4,3,11,opt{:});
imagesc( (Fsub{X}{F}(M(2)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
%print( ExampleKymographs, [saveDir,'AWAKE_CNO_30MIN.tif'], '-dtiff', '-r200' );

% ~60 mins post
subtightplot(4,3,12,opt{:});
imagesc( (Fsub{X}{F}(M(3)).cell./baseFluor)' );
set(gca,'Ytick',[], 'Xtick',[]);
xlim(frameLim);
caxis manual; caxis(colorLims);
colorbar;
%print( ExampleKymographs, [saveDir,'AWAKE_CNO_60MIN.tif'], '-dtiff', '-r200' );
print( ExampleKymographs, [saveDir,'ExampleKymographs.tif'], '-dtiff', '-r200' );
%impixelinfo;
