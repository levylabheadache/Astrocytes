clear; clc; close all;
% Load the master spreadsheet
mainDir = 'D:\2photon\'; %'D:\Andy\AstroCa\';
exptXlsx = FileFind( mainDir, 'xlsx', false, @(x)(contains( x, 'Experiment Log' )) );
exptTable = readtable( exptXlsx{1,2}, 'sheet', 'Experiments');
% Find the specific experiments needed
setGroup = {'Awake/Baseline', 'Awake/Saline', 'Awake/CNO/Ctrl', 'Awake/CLZ', 'Low CNO', 'Medium CNO', 'High CNO', 'Awake/C21', }; %{'Medium CNO', 'Medium CNO IP'}; % {'Iso/NTG', 'Awake/NTG'}; %  'Low CNO', 'Medium CNO', 'High CNO', 'Awake/NTG'}; %   {'C21'}; %  'Low CNO', 'Medium CNO', 'High CNO'    'Iso/Saline','Iso/Medium CNO','Iso/C21', 'Iso/Saline','Iso/Medium CNO', 'Saline', 'Medium CNO'
Ngroup = numel(setGroup);
for g = flip(1:Ngroup), groupRows{g} = find( strcmp( exptTable.Group, setGroup{g} ) )'; end
NgroupExpt = cellfun(@numel, groupRows);
allRows = [groupRows{:}];
Nexpt = numel( allRows );  % sum( NgroupExpt ); % size(exptTable,1);
% Initialize variables
mouse = cell(1,Nexpt); exptDate = cell(1,Nexpt); sex = cell(1,Nexpt); fov = cell(1,Nexpt); Nfov = zeros(1,Nexpt);  dataName = cell(1,Nexpt); dataDir = cell(1,Nexpt);
group = cell(1,Nexpt); drug = cell(1,Nexpt); dose = cell(1,Nexpt);
fovMovies = cell(1,Nexpt); NfovMovies = cell(1,Nexpt); LastBaseMovie = nan(1,Nexpt); FirstStimMovie = nan(1,Nexpt); LastMovie = nan(1,Nexpt);
metadata = cell(1,Nexpt); mType = cell(1,Nexpt); % mBase = cell(1,Nexpt); mStim = cell(1,Nexpt); mExpt = cell(1,Nexpt); mAcute = cell(1,Nexpt); 
baseDur = cell(1,Nexpt); baseDurTot = nan(1,Nexpt); stimDur = cell(1,Nexpt); stimDurTot = nan(1,Nexpt);
Tadj = cell(1,Nexpt); Tmin = cell(1,Nexpt); Tmed = cell(1,Nexpt); movieDur = cell(1,Nexpt);
Ncell = cell(1,Nexpt); velocity = cell(1,Nexpt); 
stimTime = cell(1,Nexpt);
proj = cell(1,Nexpt); cellROI = cell(1,Nexpt); somaROI = cell(1,Nexpt); procROI = cell(1,Nexpt); backROI = cell(1,Nexpt);
Fastro = cell(1,Nexpt); Fback = cell(1,Nexpt); dFF = cell(1,Nexpt); Fsub = cell(1,Nexpt); FsubBase = cell(1,Nexpt); C = cell(1,Nexpt); A = cell(1,Nexpt);
fSpec = cell(1,Nexpt); Pspec = cell(1,Nexpt); LFfrac = cell(1,Nexpt); Tgram = cell(1,Nexpt); fGram = cell(1,Nexpt); Pgram = cell(1,Nexpt);
cellEvt = cell(1,Nexpt); cellEvtSumm = cell(1,Nexpt); % aquaResults = cell(1,Nexpt); 
locoBout = cell(1,Nexpt); Nbout = cell(1,Nexpt); periBout = cell(1,Nexpt);  periFrame = cell(1,Nexpt); basePeri = cell(1,Nexpt); stimPeri = cell(1,Nexpt);
%bvDiam = cell(1,Nexpt);
% Gather the relevant data from each experiment
parfor x = 1:Nexpt  %1:Nexpt
    r = allRows(x);
    % Get experiment details
    mouse{x} = exptTable.Mouse{r};
    sex{x} = exptTable.Sex{r};
    exptDate{x} = num2str(exptTable.Date(r));
    group{x} = exptTable.Group{r};
    drug{x} = exptTable.Drug{r};
    dose{x} = exptTable.Dose_mg_kg_(r);    
    LastBaseMovie(x) = exptTable.LastBaseMovie(r);
    FirstStimMovie(x) = exptTable.FirstStimMovie(r);
    LastMovie(x) = min( exptTable.LastStimMovie(r), Inf );
    % Read which FOV to analyze
    tempFOV = exptTable.FOV{r}; 
    tempLims = [1, strfind( tempFOV, ', ' ), length(tempFOV)];
    for q = flip(1:numel(tempLims)-1),  fov{x}(q) = str2num(tempFOV(tempLims(q):tempLims(q+1)));  end
    Nfov(x) = numel(fov{x});

    % Load the relevant data and analyze
    for f = flip(1:Nfov(x)) %fov{X}
        % Determine which movies to analyze
        [dataName{x}{f}, dataDir{x}{f}] = SetDataPath( mouse{x}, exptDate{x}, [], fov{x}(f) );
        [metadata{x}{f},~,tempOther] = LoadProcessed( dataDir{x}{f}, 'pool' ); 
        fovMovies{x}{f} = tempOther(1).fovMovies;
        NfovMovies{x}(f) = numel(fovMovies{x}{f});
        mType{x}(f).base = find( fovMovies{x}{f} <= LastBaseMovie(x) ); % all pre-stimulus movies
        mType{x}(f).stim = find( fovMovies{x}{f} >= FirstStimMovie(x) & fovMovies{x}{f} <= LastMovie(x) ); % all relevant post-stimulus movies
        mType{x}(f).expt = [mType{x}(f).base, mType{x}(f).stim]; % all movies from the experiment (for that FOV)
        mType{x}(f).acute = [find( fovMovies{x}{f} == LastBaseMovie(x) ), find( fovMovies{x}{f} == FirstStimMovie(x) )]; % movies just before/after injection
        % total time (seconds) recorded under baseline/stimulation conditions for each FOV
        baseDur{x}(f) = sum( [metadata{x}{f}(mType{x}(f).base).duration] ); 
        stimDur{x}(f) = sum( [metadata{x}{f}(mType{x}(f).stim).duration] );
    end
    % Determine stimulus timestamp (set to end of last recording if no stimulus was applied)
    if strcmpi( exptTable.DeliveryTime{r}, 'n/a' )
        allMeta = [metadata{x}{:}]; allStamp = [allMeta.timestamp];
        [lastTime, lastTimeInd] = max( allStamp );
        stimTime{x} = lastTime + seconds( allMeta(lastTimeInd).duration ); 
    else
        stimTime{x} = datetime(exptTable.DeliveryTime{r}); 
    end
    
    % Extract and process the data
    for f = flip(1:Nfov(x)) %fov{X}
        % ASTROCYTES
        % Get astrocyte fluorescence data
        %MakeAstroROI( mouse{x}, exptDate{x}, [], fov{x}(f), 'save',true );
        [Fastro{x}{f}, Fback{x}{f}, Tadj{x}{f}, velocity{x}{f}] = GetAstroFluor(dataDir{x}{f}, 'show',false, 'refTime', stimTime{x} );  
        % Massage data: convert sec-> min, F -> dF/F, get duration and number of astrocytes
        [dFF{x}{f}, Fsub{x}{f}, FsubBase{x}{f}, C{x}{f}, A{x}{f}] = NormalizeFluor( Fastro{x}{f}, Fback{x}{f}, 'show',false );
        Ncell{x}(f) = size(Fastro{x}{f}(1).cell,2);
        Tmin{x}{f} = cellfun(@(y)(y/60), Tadj{x}{f}, 'UniformOutput',false); % convert seconds to minutes
        movieDur{x}{f} = cellfun( @range, Tmin{x}{f});
        Tmed{x}{f} = cellfun( @median, Tmin{x}{f} )';
        %PlotAstroFluor( Tmin{x}{f}, dFF{x}{f}, {}, velocity{x}{f}, 'unit','min', 'title', dataName{x}{f} ); % , Fback astroFluorFig = 
        %ReviewAstroROI(dataDir{x}{f})
        % Fourier analysis
        %{
        [fSpecTemp, PspecTemp, LFfracTemp, fGramTemp, TgramTemp, PgramTemp] = cellfun( @CaFourier, Tadj{x}{f}, {dFF{x}{f}.cell}, 'UniformOutput',false );
        fSpec{x}{f} = fSpecTemp; Pspec{x}{f} = PspecTemp; LFfrac{x}{f} = vertcat(LFfracTemp{:}); fGram{x}{f} = fGramTemp; Tgram{x}{f} = TgramTemp; Pgram{x}{f} = PgramTemp;
        %}
        % AQuA analysis
        %{
        cellEvt{x}{f} = GetAquaResults( dataDir{x}{f} ); % aquaResults{x}{f}(m)
        if ~isempty(cellEvt{x}{f})
            cellEvtSumm{x}{f}.rate = reshape([cellEvt{x}{f}.N], [NfovMovies{x}(f), Ncell{x}(f)])./repmat( movieDur{x}{f}', 1, Ncell{x}(f) ); % events/min
            cellEvtSumm{x}{f}.area = cellfun( @mean, reshape( {cellEvt{x}{f}.area}, [NfovMovies{x}(f), Ncell{x}(f)] ) ); % mean event size (square um)
            cellEvtSumm{x}{f}.dur = cellfun( @mean, reshape( {cellEvt{x}{f}.dur}, [NfovMovies{x}(f), Ncell{x}(f)] ) ); % mean event duration (seconds)
            cellEvtSumm{x}{f}.dFF = cellfun( @mean, reshape( {cellEvt{x}{f}.dFF}, [NfovMovies{x}(f), Ncell{x}(f)] ) ); % mean event peak dF/F      
        else
            fprintf('\nNo AQuA data for %s\n', dataName{x}{f} );
        end
        %}
        
        % LOCOMOTION
        % {
        for m = 1:NfovMovies{x}(f) % mExpt{x}{f}
            [locoBout{x}{f}{m}, Nbout{x}{f}(m)] = LocoBouts(Tadj{x}{f}{m}, velocity{x}{f}{m}, 'show',false, 'pause',false ); % 
            [periBout{x}{f}(m), periFrame{x}{f}{m}] = GetPeriBout( Tadj{x}{f}{m}, velocity{x}{f}{m}, locoBout{x}{f}{m}, dFF{x}{f}(m), A{x}{f}(m), 'base',10, 'run',10, 'inspect',false, 'show',false ); %  Fsub{x}{f}(m)
        end
        %}
        % VASCULATURE
        %{
        %MakeVascROI( mouse{x}, exptDate{x}, [], fov{x}(f), 'save',true );
        %[bvDiam{x}{f}, ~] = GetVascDiameter(dataDir{x}{f}, 'show',true, 'refTime', exptTable.DeliveryTime(r));
        %}
    end
    baseDurTot(x) = sum( baseDur{x} ); % total baseline recording time (in seconds) over all FOV
    stimDurTot(x) = sum( stimDur{x} );
end
Nastro = cellfun( @sum, Ncell );
xGroup = cell(1,Ngroup);
for g = 1:Ngroup,  xGroup{g} = find( strcmp( group, setGroup{g} ) );  end 
groupN = cellfun( @numel, xGroup );
% Combine all data, across fields of view, from each experiment an into a set of supermatrices
[Tcomb, fovComb, Fcomb, dFFcomb, Acomb, locoComb, fovCol] = cellfun( @CombineFOVdata, Tmin, Fsub, dFF, A, Ncell, mType, velocity, locoBout, periFrame, 'UniformOutput',false ); %CombineFOVdata( Tmin{x}, Fastro{x}, Ncell{x}, mExpt{x} )
[Tmid, binSumm, Nbin, binInt] = cellfun( @BinCombinedData, Tcomb, fovComb, Fcomb, dFFcomb, Acomb, locoComb, 'UniformOutput',false ); binInt = binInt{1};

%BinCombinedData(Tcomb{x}, fovComb{x}, Fcomb{x}, dFFcomb{x}, Acomb{x}, locoComb{x}, 'show',true, 'interval',0.5);

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
