%% Pool summary statistics by group
FmedCell = cell(1,Nexpt); FmedSoma = cell(1,Nexpt); FmedProc = cell(1,Nexpt);
AbinCell = cell(1,Nexpt); AbinSoma = cell(1,Nexpt); AbinProc = cell(1,Nexpt);
LFratio = cell(1,Nexpt);
travel = cell(1,Nexpt); runTime = cell(1,Nexpt);
for x = 1:Nexpt %xGroup{g}
    Fmed = [binSumm{x}.Fcen]; %  [binSumm{x}.Fmed]; %  
    FmedCell{x} = vertcat(Fmed.cell);
    FmedSoma{x} = vertcat(Fmed.soma);
    FmedProc{x} = vertcat(Fmed.proc);
    
    Abin = [binSumm{x}.Abin]; % 
    AbinCell{x} = vertcat(Abin.cell);
    AbinSoma{x} = vertcat(Abin.soma);
    AbinProc{x} = vertcat(Abin.proc);
    
    LFratio{x} = vertcat( binSumm{x}.LFratio );
    travel{x} = vertcat( binSumm{x}.travel );
    runTime{x} = vertcat( binSumm{x}.runTime );
end

CountColNumbers = @(x)( sum(~isnan(x),1) );
CountRowNumbers = @(x)( sum(~isnan(x),2) );
baseBin = find( Tmid{1} < 0 );
groupFcell = cell(1,Ngroup); groupFsoma = cell(1,Ngroup); groupFproc = cell(1,Ngroup); 
groupAcell = cell(1,Ngroup); groupAsoma = cell(1,Ngroup); groupAproc = cell(1,Ngroup); 
groupFnorm = cell(1,Ngroup); groupNastro = zeros(1,Ngroup);
groupLF = cell(1,Ngroup); groupSP = cell(1,Ngroup);  
groupTravel = cell(1,Ngroup); groupRun = cell(1,Ngroup);
for g = 1:Ngroup
    groupFcell{g} = horzcat(FmedCell{xGroup{g}});
    groupFnorm{g} = groupFcell{g}./nanmean(groupFcell{g}(baseBin,:), 1);
    groupFsoma{g} = horzcat(FmedSoma{xGroup{g}});
    groupFproc{g} = horzcat(FmedProc{xGroup{g}});
    groupAcell{g} = horzcat(AbinCell{xGroup{g}});
    groupAsoma{g} = horzcat(AbinSoma{xGroup{g}});
    groupAproc{g} = horzcat(AbinProc{xGroup{g}});
    groupSP{g} = groupFsoma{g}./groupFproc{g};
    groupNastro(g) = sum(Nastro(xGroup{g}));
    groupLF{g} = horzcat( LFratio{xGroup{g}} );
    groupTravel{g} = horzcat( travel{xGroup{g}} );
    groupRun{g} = horzcat( runTime{xGroup{g}} );
end
groupFnormMean = cellfun(@nanmean, groupFnorm, repmat({2}, 1, Ngroup), 'UniformOutput',false);
groupFnormStd = cellfun(@nanstd, groupFnorm, repmat({[]}, 1, Ngroup), repmat({2}, 1, Ngroup), 'UniformOutput',false);
groupBinNastro = cellfun(CountRowNumbers, groupFnorm, 'UniformOutput',false);
groupBinNmouse = cellfun(CountRowNumbers, groupTravel, 'UniformOutput',false);

%% Visualize the times at which each group was most-observed
close all; figure('WindowState','max');
for g = 1:Ngroup
    subplot(Ngroup,1,g);
    plot( Tmid{1}, groupBinNastro{g}' ); hold on;
    ylabel('Cells Observed');
    title( setGroup{g} );
end
xlabel('Time Bin');

%%
for x = 1:Nexpt
    for f = 1:Nfov(x)
        PlotAstroFluor( Tmin{x}{f}, dFF{x}{f}, {}, velocity{x}{f}, 'unit','min', 'title', dataName{x}{f} );
        pause;
        close all;
    end
end


%% Show effect on NORMALIZED FLUOR, each experiment within each group
close all;
figure('WindowState','max');
for g = 1:Ngroup
    c = 0; clearvars sp; clf;
    for x = xGroup{g}
        fovFnorm = FmedCell{x}./nanmean(FmedCell{x}(baseBin,:), 1);
        tempNastro = CountRowNumbers(fovFnorm);
        tempMean = nanmean(fovFnorm,2);
        tempSEM = nanstd(fovFnorm,0,2)./sqrt(tempNastro);
        c = c+1;
        sp(c) = subplot(NgroupExpt(g), 1, c);
        %tempNastroFovStr = ;
        line([-60,120], [1,1], 'Color','r', 'LineStyle','--'); hold on;
        plot( Tmid{1}, fovFnorm, '.' ); hold on;
        errorbar( Tmid{1},  tempMean, tempSEM, 'CapSize',3, 'Color','k' ); hold on;
        title( sprintf('%s [g,x] = [%i, %i]:  %s %s (%i FOV, [%s] astrocytes respectively)',setGroup{g}, g, x, mouse{x}, exptDate{x}, Nfov(x), sprintf('%i, ', Ncell{x} )) );
        set(gca,'Xtick',-60:10:120, 'TickDir','out', 'TickLength',[0.005,0]);
        ylabel('Norm. Fluor');
    end
    linkaxes(sp,'xy');
    if g < Ngroup, pause; end
end


%% Show effect on RUN TIME, each experiment within each group
close all;
figure('WindowState','max');
for g = 1:Ngroup
    c = 0; clearvars sp; clf;
    for x = xGroup{g}
        c = c+1;
        sp(c) = subplot(NgroupExpt(g), 1, c);
        plot( Tmid{1}, runTime{x}, '.' ); hold on;
        title( sprintf('%s [g,x] = [%i, %i]:  %s %s (%i FOV)',setGroup{g}, g, x, mouse{x}, exptDate{x}, Nfov(x)) );
        set(gca,'Xtick',-60:10:120, 'TickDir','out', 'TickLength',[0.005,0]);
        ylabel('Run Time (s)');
        ylim([0, interval*60]);
    end
    linkaxes(sp,'xy');
   	if g < Ngroup, pause; end
end


%% Show overall effect on NORMALIZED FLUOR, by group
opt = {[0.08,0.04], [0.08,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
figure('WindowState','max');
for g = 1:Ngroup
    subtightplot(Ngroup,1,g,opt{:});
    line([-60,120], [1,1], 'Color','r', 'LineStyle','--'); hold on;
    errorbar( Tmid{1}, groupFnormMean{g}, groupFnormStd{g}./sqrt(groupBinNastro{g}), 'CapSize',3 ); hold on;
    %set(gca,'CapSize',5);
    ylim([0.8,2]);
    title( sprintf('%s (%i experiments, %i astrocytes)', setGroup{g}, NgroupExpt(g), groupNastro(g) ), 'FontSize',10);
    ylabel('Norm. Med. Fcell');
    %pause;
end

%% Plot soma and process fluor, by group
close all; clearvars sp;
figure('WindowState','max');
opt = {[0.08,0.04], [0.08,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
for g = flip(1:Ngroup)
    sp(g) = subtightplot(Ngroup,1,g,opt{:});
    %yyaxis left
    errorbar( Tmid{1}, nanmean(groupFsoma{g},2), nanstd(groupFsoma{g},[],2)./sqrt(groupBinNastro{g}) ); hold on;
    errorbar( Tmid{1}, nanmean(groupFproc{g},2), nanstd(groupFproc{g},[],2)./sqrt(groupBinNastro{g}) );
    title( sprintf('%s (%i experiments, %i astrocytes, %2.1f minute intervals)', setGroup{g}, NgroupExpt(g), groupNastro(g), interval ), 'FontSize',12);
    xlim([-Inf,Inf]);  %ylim([1, 2]);
    ylabel('Med. Fluor.'); 
    if g == 1, legend('Soma','Processes'); end
    set(gca, 'FontSize',12 )
    %pause;
end
xlabel('Time (min)');
linkaxes(sp,'x');


%% Plot activity, by group
close all; clearvars sp;
figure('WindowState','max');
opt = {[0.08,0.04], [0.08,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
for g = flip(1:Ngroup)
    sp(g) = subtightplot(Ngroup,1,g,opt{:});
    groupMean{g} = nanmean(groupAcell{g},2); % groupSEM{g}(isnan(groupSEM{g})) = 0;
    groupSEM{g} = nanstd(groupAcell{g},[],2)./sqrt(groupBinNastro{g}); %groupSEM{g}(isnan(groupSEM{g})) = 0;
    %plot( Tmid{1}, groupMean{g}, '.' );
    errorbar( Tmid{1}, groupMean{g}, groupSEM{g} ); hold on;
    title( sprintf('%s (%i experiments, %i astrocytes, %2.1f minute intervals)', setGroup{g}, NgroupExpt(g), groupNastro(g), binInt ), 'FontSize',12);
    xlim([-Inf,Inf]);  %ylim([1, 2]);
    ylabel('Med. Act.'); 
    if g == 1, legend('Soma','Processes'); end
    set(gca, 'FontSize',12 )
    %pause;
end
xlabel('Time (min)');
linkaxes(sp,'x');

%% Plot soma-process fluor ratio, by group
close all; clearvars sp;
figure('WindowState','max');
opt = {[0.06,0.04], [0.08,0.04], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
for g = 1:Ngroup
    sp(g) = subtightplot(Ngroup,1,g,opt{:});
    errorbar( Tmid{1}, nanmean(groupSP{g},2), nanstd(groupSP{g},[],2)./sqrt(groupBinNastro{g}) ); hold on;
    title( sprintf('%s (%i experiments, %i astrocytes, %2.1f minute intervals)', setGroup{g}, NgroupExpt(g), groupNastro(g), interval ), 'FontSize',12);
    %xlim([-Inf,Inf]);  %ylim([1, 2]);
    ylabel('Med. Fluor.'); 
    set(gca, 'FontSize',12 )
    %pause;
end
xlabel('Time (min)');
linkaxes(sp,'xy');

%% Plot low-frequency power ratio by group
close all; clearvars sp;
figure('WindowState','max');
for g = flip(1:Ngroup)
    sp(g) = subplot(Ngroup,1,g); 
    errorbar( Tmid{1}, nanmean(groupLF{g},2), nanstd(groupLF{g},[],2)./sqrt(groupBinNastro{g}) ); hold on;
    title( sprintf('%s (%i experiments, %i astrocytes, %2.1f minute intervals)', setGroup{g}, NgroupExpt(g), groupNastro(g), interval ), 'FontSize',10);
    ylabel('LF Power Ratio');
    ylim([0,1]);
end
linkaxes(sp,'xy');

%% Plot wheel travel, by group
close all; clearvars sp;
figure('WindowState','max');
for g = 1:Ngroup
    sp(g) = subplot(Ngroup,1,g); 
    errorbar( Tmid{1}, nanmean(groupTravel{g},2), nanstd(groupTravel{g},[],2)./sqrt(groupBinNmouse{g}) ); hold on;
    title( sprintf('%s (%i experiments, %2.1f minute intervals)', setGroup{g}, NgroupExpt(g), interval ), 'FontSize',10);
    ylabel('Distance Traveled (cm)');
    %ylim([0,1]);
end
linkaxes(sp,'xy');

%% Plot wheel run time, by group
close all; clearvars sp;
figure('WindowState','max');
for g = 1:Ngroup
    sp(g) = subplot(Ngroup,1,g); 
    errorbar( Tmid{1}, nanmean(groupRun{g},2), nanstd(groupRun{g},[],2)./sqrt(groupBinNmouse{g}) ); hold on;
    title( sprintf('%s (%i experiments, %2.1f minute intervals)', setGroup{g}, NgroupExpt(g), interval ), 'FontSize',10);
    ylabel('Run Time(s)');
    %ylim([0,1]);
end
linkaxes(sp,'xy');

%% Compare specific pre-injection, acute and sustained timepoints
preMean = cell(1,Ngroup); preStd = cell(1,Ngroup); NpreObs = cell(1,Ngroup);
acuteMean = cell(1,Ngroup); acuteStd = cell(1,Ngroup); NacuteObs = cell(1,Ngroup);
sustMean = cell(1,Ngroup); sustStd = cell(1,Ngroup); NsustObs = cell(1,Ngroup);

preBins = find( Tmid{1} >= -20 & Tmid{1} <= 0 );
acuteBins = find( Tmid{1} >= 10 & Tmid{1} <= 30 );
sustBins = find( Tmid{1} >= 50 & Tmid{1} <= 70 );
for g = 1:Ngroup
    preFcell = groupFcell{g}(preBins,:);  
    NpreObs{g} = CountColNumbers( preFcell ); % how many observations were made during pre-injection window
    preMean{g} = mean( preFcell, 1, 'omitnan' );
    %preStd{g} = std( preFcell, 0, 1, 'omitnan' );
    preMean{g}(NpreObs{g} < 3) = NaN;
    
    acuteFcell = groupFcell{g}(acuteBins,:);
    NacuteObs{g} = CountColNumbers( acuteFcell ); % how many observations were made during acute window
    acuteMean{g} = mean( acuteFcell, 1, 'omitnan' );
    %acuteStd{g} = std( acuteFcell, 0, 1, 'omitnan' );
    acuteMean{g}(NacuteObs{g} < 3) = NaN;

    sustFcell = groupFcell{g}(sustBins,:);
    sustMean{g} = mean( sustFcell, 1, 'omitnan' );
    %sustStd{g} = std( sustFcell, 0, 1, 'omitnan' );
    NsustObs{g} = CountColNumbers( sustFcell ); % how many observations were made during sustained window
    sustMean{g}(NsustObs{g} < 3) = NaN;
    
    pasMat{g} = vertcat( preMean{g}, acuteMean{g}, sustMean{g} );
    pasMat{g} = pasMat{g}(:,~isnan( sum( pasMat{g}, 1 ) ));   % Remove incomplete columns
end

% Unnormalized
close all; clearvars sp;
figure('WindowState','max');
for g = 1:Ngroup
    sp(g) = subplot(1,Ngroup,g);
    plot( pasMat{g} ); hold on;
    errorbar( 1:3, nanmean( pasMat{g}, 2 ), nanstd( pasMat{g}, 0, 2)./sqrt(CountRowNumbers(pasMat{g})), '-.', 'color','k','MarkerSize',6, 'LineWidth',3 )

    xlim([0.5,3.5]);
    set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Sust'});
    title( sprintf('%s', setGroup{g}) );
end

% Normalized to pre-injection
close all; clearvars sp;
figure('WindowState','max');
for g = 1:Ngroup
    pasNorm{g} = (pasMat{g}-pasMat{g}(1,:))./pasMat{g}(1,:);
    sp(g) = subplot(1,Ngroup,g);
    plot( pasNorm{g} ); hold on;
    errorbar( 1:3, nanmean( pasNorm{g}, 2 ), nanstd(pasNorm{g},0,2)./sqrt(CountRowNumbers( pasNorm{g} )), '-.', 'color','k','MarkerSize',6, 'LineWidth',3 ); % SEM(pasNorm{g})

    xlim([0.5,3.5]);
    set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Sust'});
    title( sprintf('%s', setGroup{g}) );
end

%% Repeated measures ANOVA on the pre/acute/sustained fluor data
pasTbl = cell(1,Ngroup);
for g = 1:Ngroup
    tempMat = pasMat{g}';
    
    pasTbl{g} = table( repmat(setGroup(g), size(tempMat,1), 1), tempMat(:,1), tempMat(:,2), tempMat(:,3), 'VariableNames', {'Group','Pre','Acute','Sustained'}); % , 'VariableNames', {'Cell','Pre','Acute','Sustained'}
    %Time = [0,30,60];
    %tempMod = fitrm(tempTbl, 't0-t2~1', 'WithinDesign',Time);
    %tempAnovaRes = ranova( tempMod );
end
fullTbl = vertcat( pasTbl{:} );
writetable( fullTbl, [mainDir,'pasData.xls'])

%% Determine how long it took for each cell to become activated and then deactivated after injection
Tact = cell(1,Ngroup); Tdeact = cell(1,Ngroup);
Nsigma = 2.5;
preInd = find( Tmid{1} > 5, 1 );
for g = 2 %1:Ngroup    
    groupFmean = nanmean(groupFcell{g}(baseBin,:), 1);
    groupFstd = nanstd(groupFcell{g}(baseBin,:), 0, 1);
    groupFthresh = groupFmean + Nsigma*groupFstd;
    Tact{g} = nan(1, groupNastro(g)); Tdeact{g} = nan(1, groupNastro(g));
    for q = 1:groupNastro(g)
        % Activation
        tempAct = find( groupFcell{g}(:,q) > groupFthresh(q) );
        if ~isempty( tempAct)
            tempAct = tempAct(tempAct > preInd );
            actInd = tempAct(1);
            Tact{g}(q) = Tmid{1}(actInd);
            % Deactivation
            tempDeact = find( groupFcell{g}(:,q) < groupFthresh(q) );
            tempDeact = tempDeact(tempDeact > actInd );
            if ~isempty(tempDeact)
                deactInd = tempDeact(1);
                Tdeact{g}(q) = Tmid{1}(deactInd);
            end
        end

        % Plot
        plot( Tmid{1}, groupFcell{g}(:,q) ); hold on;
        line([-60,120], groupFmean(q)*[1,1], 'color','k');
        line([-60,120], groupFthresh(q)*[1,1], 'color','r','linestyle','--');
        title( sprintf('%2.1f, %2.1f', Tact{g}(q), Tdeact{g}(q)) );
        pause;
        hold off;
    end
end

%% Show effect on ACTIVATION, each experiment within each group
actFrac = cell(1,Ngroup); %Tdeact = cell(1,Ngroup);
Tact = cell(1,Ngroup); Tdeact = cell(1,Ngroup);
Nsigma = 3;
preInd = find( Tmid{1} > 5, 1 );
postBin = find( Tmid{1} > 5 );
close all;
figure('WindowState','max');
for g = 1:Ngroup
    Tact{g} = nan(1,NgroupExpt(g)); Tdeact{g} = nan(1,NgroupExpt(g));
    c = 0; clearvars sp; clf;
    for x = xGroup{g}
        c = c+1;
        fovFnorm = FmedSoma{x}./nanmean(FmedSoma{x}(baseBin,:), 1); % FmedProc{x}./nanmean(FmedProc{x}(baseBin,:), 1);   % FmedCell{x}./nanmean(FmedCell{x}(baseBin,:), 1);
        tempNastro = CountRowNumbers(fovFnorm);
        tempMean = nanmean(fovFnorm,2);
        baseMean = nanmean( tempMean(baseBin));
        baseStd = nanstd(tempMean(baseBin));
        tempThresh = baseMean + Nsigma*baseStd;
        postMean = tempMean(postBin);
        
        % Activation - fraction of time spent above threshold
        Nobs = sum(~isnan(postMean));
        Nact = sum(postMean > tempThresh);
        actFrac{g}(c) = Nact/Nobs;
        % Activation - using threshold crossings only
        %{
        tempAct = find( tempMean > tempThresh );
        tempAct = tempAct(tempAct > preInd );
        if ~isempty( tempAct)
            actInd = tempAct(1);
            Tact{g}(c) = Tmid{1}(actInd);
            % Deactivation
            tempDeact = find( tempMean < tempThresh );
            tempDeact = tempDeact(tempDeact > actInd );
            if ~isempty(tempDeact)
                deactInd = tempDeact(1);
                Tdeact{g}(c) = Tmid{1}(deactInd);
            end
        end
        %}
        % Plot
        sp(c) = subplot(NgroupExpt(g), 1, c);
        line([-60,120], baseMean*[1,1], 'Color','k'); hold on;
        line([-60,120], tempThresh*[1,1], 'Color','r', 'LineStyle','--'); hold on;
        plot( Tmid{1}, tempMean ); hold on;
        title( sprintf('%s [g,x] = [%i, %i]:  %s %s (%i FOV, [%s] astrocytes respectively): Activation fraction = %2.2f',setGroup{g}, g, x, mouse{x}, exptDate{x}, Nfov(x), sprintf('%i, ',Ncell{x}), actFrac{g}(c) ) );
        set(gca,'Xtick',-60:10:120, 'TickDir','out', 'TickLength',[0.005,0]);
        ylabel('Norm. Fluor');
        %title( sprintf('%2.1f, %2.1f', Tact{g}(c), Tdeact{g}(c)) );
    end
    linkaxes(sp,'xy');
    if g < Ngroup, pause; end
end
%% Show effect on ACTIVATION, each ASTROCYTE FROM EACH FOV within each group
Nsigma = 2.5;
onWindow = find( Tmid{1} > 5 & Tmid{1} < 20 );
linThresh = 30; [~,linLim] = min( abs(Tmid{1} - linThresh) );
longThresh = 90; 
%[~,offLim] = min( abs(Tmid{1} - 90) );
close all; figure('WindowState','max');
Ton = cell(1,Nexpt); Toff = cell(1,Nexpt); longAct = cell(1,Nexpt); linSlope = cell(1,Nexpt);
for g = 4 % 1:Ngroup
    for x = xGroup{g}
        Ton{x} = cell(1,Nfov(x)); Toff{x} = cell(1,Nfov(x));
        for f = find( lastObs{x} >= 60 ) % 1:Nfov(x)
            Ton{x}{f} = nan(1,Ncell{x}(f));  Toff{x}{f} = nan(1,Ncell{x}(f)); longAct{x}{f} = nan(1,Ncell{x}(f)); linSlope{x}{f} = nan(1,Ncell{x}(f));
            fovFnorm = FmedCell{x}(:,fovCol{x}{f})./nanmean(FmedCell{x}(baseBin,fovCol{x}{f}), 1); 
            baseMean = nanmean( fovFnorm(baseBin,:), 1);
            baseStd = nanstd( fovFnorm(baseBin,:),0, 1 );
            actThresh = baseMean+Nsigma*baseStd;
            % For each cell, find the timepoints when the fluorescence crosses the threshold
            for a = 1:Ncell{x}(f)
                tempActBin = double( fovFnorm(:,a) > actThresh(a) )';
                tempActBin( isnan(fovFnorm(:,a)) ) = NaN;
                onBin = strfind( tempActBin, [1,1] ); 
                if ~isempty( onBin ), onBin = onBin( ismember( onBin, onWindow ) ); end
                if ~isempty( onBin )
                    onBin = onBin(1);
                else
                    onBin = NaN; offBin = NaN;
                end
                if ~isnan(onBin)
                    offBin = strfind( tempActBin, [0,0] ); 
                    offBin = offBin( offBin > onBin );
                    if ~isempty( offBin )
                        offBin = offBin(1);
                    else
                        offBin = find( tempActBin == 1, 1, 'last' ); % 
                        %offBin(offBin > offLim ) = []; % enforce a maximum to compensate for inconsistent duration of observation
                        %offBin = offBin(end);
                    end
                    Ton{x}{f}(a) = Tmid{1}(onBin);
                    Toff{x}{f}(a) = Tmid{1}(offBin);
                    if Toff{x}{f}(a) >= longThresh, longAct{x}{f}(a) = 1;  else, longAct{x}{f}(a) = 0; end
                end
                % Linear regression of the post-acute phase
                linObsBin = find( ~isnan(fovFnorm(:,a)) );
                linObsBin = linObsBin(linObsBin >= linLim );
                if numel(linObsBin > 5)
                    Tlin = Tmid{1}( linObsBin )';
                    Flin = fovFnorm( linObsBin, a );
                    p = polyfit( Tlin, Flin, 1 );
                    linSlope{x}{f}(a) = p(1);
                end
                %Tlin = 
                
                
                % {
                clf;
                line([-60,120], actThresh(a)*[1,1], 'Color','k', 'LineStyle','--'); hold on;
                plot( Tmid{1}, fovFnorm(:,a) ); hold on;
                plot( Tlin, Flin, 'ko' )
                if ~isnan(onBin),  plot( Ton{x}{f}(a), fovFnorm(onBin,a), 'ko' );  end
                if ~isnan(offBin),  plot( Toff{x}{f}(a), fovFnorm(offBin,a), 'kx' ); end
                title( sprintf('%s [g,x,f,a] = [%i, %i, %i, %i]:  %s %s FOV %i - Ton = %2.1f,  Toff = %2.1f, longAct = %1.1f', setGroup{g}, g, x, f, a,  mouse{x}, exptDate{x}, fov{x}(f), Ton{x}{f}(a), Toff{x}{f}(a), longAct{x}{f}(a)) );
                set(gca,'Xtick',-60:10:120, 'TickDir','out', 'TickLength',[0.005,0]);
                ylabel('Norm. Fluor');
                %title( sprintf('%2.1f, %2.1f', Tact{g}(c), Tdeact{g}(c)) );
                xlim([-60,120]);
                %hold off;
                pause; %clf;
                %}
            end
        end
    end
end



%% Compare FRACTION OF CELLS ACTIVATED, AND DEMONSTRATING LONG-LASTING ACTIVATION

for g = 1:Ngroup
    tempOn = [Ton{xGroup{g}}]; tempOn = [tempOn{:}];
    tempLong = [longAct{xGroup{g}}]; tempLong = [tempLong{:}];
    Non = CountRowNumbers( tempOn ); % how many cells became activated?
    Nlong = sum( tempLong == 1 );
    onFrac(g) = Non/length(tempOn);
    longFrac(g) = Nlong/Non;
end

figure('WindowState','maximized');
subplot(1,2,1);
bar( onFrac )
ylabel('Fraction of Cells Activated');
xlim([0, Ngroup]+0.5);
set(gca,'Xtick',1:Ngroup, 'XtickLabel',setGroup, 'box','off');
xtickangle(30); axis square;
subplot(1,2,2);
bar( longFrac )
ylabel('Fraction of Activated Cells');
xlim([0, Ngroup]+0.5);
set(gca,'Xtick',1:Ngroup, 'XtickLabel',setGroup, 'box','off');
xtickangle(30); axis square; 

%% Compare ACTIVATION DURATION
figure;
groupActDur = cell(1,Ngroup);
for g = 1:Ngroup
    tempOn = [Ton{xGroup{g}}]; tempOn = [tempOn{:}];
    tempOff = [Toff{xGroup{g}}]; tempOff = [tempOff{:}];
    groupDur{g} = tempOff - tempOn;
    plot( g, groupDur{g}, '.', 'MarkerSize',10 ); hold on;
    errorbar( g, nanmean(  groupDur{g} ), SEM( groupDur{g} ) );
end
xlim([0, Ngroup]+0.5);

%% Example Kymographs
Tpre = [-20,0]; Tacute = [20,40]; Tsust = [50,70];

preBins = find( Tmid{1} >= Tpre(1) & Tmid{1} <= Tpre(2) );
acuteBins = find( Tmid{1} >= Tacute(1) & Tmid{1} <= Tacute(2) );
sustBins = find( Tmid{1} >= Tsust(1) & Tmid{1} <= Tsust(2) );
minObs = 2;
chunkLength = numel(preBins);
chunk = cell(Ngroup,1); Ngood = nan(1,Ngroup); %nan(Ngroup,3); % 
figure;
for g = 1:Ngroup
    preChunk = groupFcell{g}(preBins,:);
    preObs = CountColNumbers( preChunk ); % how many bins were observed for each cell?
    acuteChunk = groupFcell{g}(acuteBins,:);
    acuteObs = CountColNumbers( acuteChunk ); % how many bins were observed for each cell?
    sustChunk = groupFcell{g}(sustBins,:);
    sustObs = CountColNumbers( sustChunk ); % how many bins were observed for each cell?
    
    % Determine well-observed cells
    goodCells = find( preObs >= minObs & acuteObs >= minObs & sustObs >= minObs );
    Ngood(g) = numel(goodCells);
    preMean = nanmean( preChunk(:,goodCells) );
    chunk{g} = vertcat( nanmean( preChunk(:,goodCells) ), nanmean( acuteChunk(:,goodCells) ), nanmean( sustChunk(:,goodCells) ) )./preMean;
    
    %imagesc( chunk{g}' );
    %{
    % Gather pre/30 min/60 min chunks
    chunk{g} = groupFcell{g}( [preBins, acuteBins, sustBins], : );
    %CountColNumbers( chunk{g} )
    
    % Remove nan columns
    % Plot the results
    imagesc( chunk{g}' );
    set(gca, 'Xtick',[1,chunkLength,chunkLength+1,2*chunkLength,2*chunkLength+1,3*2*chunkLength], 'XtickLabels', horzcat( Tpre, Tacute, Tsust ), 'TickDir','out', 'box','off' );
    pause;
    %}
end
groupEdge = [0, cumsum(Ngood)] + 0.5;
allChunk = horzcat( chunk{:})';
allLength = size(allChunk,1);

% Make the figure
FS = 14;
figure('WindowState','maximized');
imagesc( allChunk ); hold on;
set(gca,'Xtick',1:3, 'XtickLabel',{'-20 - 0', '20 - 40', '50 - 70'}, 'Ytick',groupEdge(1:end-1) + Ngood/2, 'YtickLabel', setGroup, 'TickDir','out', 'FontSize',FS );
cb = colorbar; 
cb.Label.String = 'Baseline-Normalized Mean Fluorescence';
cb.FontSize = FS;
line([1.5,1.5],[1,allLength], 'Color','k' ); 
line([2.5,2.5],[1,allLength], 'Color','k' );
for g = 1:Ngroup
    line([0.5,3.5], groupEdge(g)*[1,1], 'Color','k' );
end
xlabel('Time Relative to Injection (mins)');
