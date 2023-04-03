function [periBout, periInd] = GetPeriBout( T, velocity, locoBout, F, A, varargin )
% PeriBoutData gather relevant data from a defined period of time around the onset of each bout of running. For use with PeriBoutStats.  
IP = inputParser;
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'velocity', @isnumeric )
addRequired( IP, 'locoBout', @isstruct )
addRequired( IP, 'F', @isstruct )
addRequired( IP, 'A', @isstruct )
addParameter( IP, 'base', 10, @isnumeric )
addParameter( IP, 'run', 15, @isnumeric )
addParameter( IP, 'init', 1, @isnumeric )
addParameter( IP, 'latLim', [0.3,5], @isnumeric )
addParameter( IP, 'testLim', [-3,3], @isnumeric )
addParameter( IP, 'Nboot', 500, @isnumeric )
addParameter( IP, 'zMin', 3, @isnumeric )
addParameter( IP, 'Ntrig', 3, @isnumeric ) % minimum # of frames, within latLim, with fluor above zMin for an event to count as run-triggered
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'inspect', false, @islogical )
parse( IP, T, velocity, locoBout, F, A, varargin{:} ); 
baseDur = IP.Results.base;
runDur = IP.Results.run;
initDur = IP.Results.init; 
testLim = IP.Results.testLim;
Nboot = IP.Results.Nboot;
show = IP.Results.show; 
inspectBout = IP.Results.inspect;

Nframe = size(T,1);
Ncell = size(F.cell, 2);
Nbout = numel(locoBout );
frameRate = 1/median(diff(T));
%baseDur = 3; runDur = 2; binDur = 0.5; minIso = 2*baseDur;
NbaseFrame = ceil(frameRate*baseDur);
NrunFrame = ceil(frameRate*runDur);
NtotFrame = NbaseFrame + NrunFrame;
NinitFrame = ceil(frameRate*initDur);
cspStruct = struct('F',nan(NtotFrame,Ncell,Nbout), 'A',nan(NtotFrame,Ncell,Nbout), 'NAD',nan(Nbout,Ncell), 'NADci', nan(Nbout,Ncell,2), ...
    'trig',false(Nbout,Ncell), 'trigFrac',nan(Nbout,1), 'thresh',nan(Nbout,Ncell), 'latency',nan(Nbout,Ncell), 'peak',nan(Nbout,Ncell), 'delta',nan(Nbout,Ncell) ); 

periBout = struct('Nbout',Nbout, 'Ncell',Ncell, 'Nframe',NtotFrame, 'base', 1:NbaseFrame, 'Nbase', NbaseFrame, 'run',NbaseFrame+1:NtotFrame, 'Nrun',NrunFrame, 'Ninit',NinitFrame,...
        'frames',nan(NtotFrame,Nbout), 'T',nan(NtotFrame,1), 'velocity',nan(NtotFrame,Nbout), 'dur',nan(Nbout,1), 'maxSpd',nan(Nbout,1), 'travel',nan(Nbout,1), 'initMax', nan(Nbout,1),...
        'cell', cspStruct', 'soma', cspStruct, 'proc',cspStruct );
periInd = struct('base',[], 'run',[], 'all',[]);
if Nbout > 0 && Nframe >= NtotFrame
    periBout.T = double(T(1:NtotFrame) - T(NbaseFrame+1));
    if inspectBout
        FS = 12;
        close all;
        figure('WindowState','max', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        sp(4) = subplot(4,1,4); sp(3) = subplot(4,1,3); sp(2) = subplot(4,1,2); sp(1) = subplot(4,1,1);
    end   
    ZtestBase = find( periBout.T < 0 & periBout.T >= testLim(1) ); %  & ~isnan(periBout.proc(:,1,b))
    ZtestRun = find( periBout.T >= 0 & periBout.T <= testLim(2) );  %  & ~isnan(periBout.proc(:,1,b))
    for b = 1:Nbout
        tempFrames = locoBout(b).frames(1)-NbaseFrame:locoBout(b).frames(1)-NbaseFrame+NtotFrame-1;
        goodZ = find(tempFrames>0 & tempFrames<Nframe+1);
        goodBase = goodZ( goodZ > 0 & goodZ <= NbaseFrame );
        goodRun = goodZ( goodZ > NbaseFrame & goodZ <= NtotFrame );
        tempFrames( tempFrames < 1 | tempFrames > Nframe ) = [];        
        tempInit = locoBout(b).frames(1):locoBout(b).frames(1)+NinitFrame-1;  tempInit( tempInit < 1 | tempInit > Nframe ) = [];
        % Define test intervals just before/after run onset
        goodTestBase = intersect( ZtestBase, goodBase ); NtestBase = numel( goodTestBase ); 
        goodTestRun = intersect( ZtestRun, goodRun ); NtestRun = numel( goodTestRun );
        % GET LOCOMOTION, FLUOR AND ACTIVITY DATA
        periBout.frames(goodZ,b) = tempFrames;
        periBout.velocity(goodZ,b) = velocity(tempFrames);
        %periBout.accel(:,b) = accel(periBout.frames(:,b));
        periBout.dur(b) = locoBout(b).dur;
        periBout.maxSpd(b) = locoBout(b).maxSpd;
        periBout.travel(b) = locoBout(b).travel;
        periBout.cell.F(goodZ,:,b) = F.cell(tempFrames,:); %
        periBout.soma.F(goodZ,:,b) = F.soma(tempFrames,:);
        periBout.proc.F(goodZ,:,b) = F.proc(tempFrames,:);
        periBout.cell.A(goodZ,:,b) = A.cell(tempFrames,:); %
        periBout.soma.A(goodZ,:,b) = A.soma(tempFrames,:);
        periBout.proc.A(goodZ,:,b) = A.proc(tempFrames,:);
        % QUANTIFY INITIAL LOCOMOTION PROPERTIES
        periBout.initMax(b) = max( abs(velocity(tempInit)) ); %max(periBout.velocity(temp,b));
        % DETERMINE WHETHER THE BOUT APPEARED TO TRIGGER CA ACTIVITY, AND QUANTIFY RESPONSE
        minNtest = round(2*min(periBout.Nbase, NtestRun)/3);
        if NtestBase >= minNtest && NtestRun >= minNtest
            % Compare activity just before/after onset of locomotion
            CellBaseTestData = periBout.cell.A(ZtestBase,:,b);
            CellRunTestData = periBout.cell.A(ZtestRun,:,b);
            periBout.cell.NAD(b,:) = SetNormDiff( CellBaseTestData, CellRunTestData );
            SomaBaseTestData = periBout.soma.A(ZtestBase,:,b);
            SomaRunTestData = periBout.soma.A(ZtestRun,:,b);
            periBout.soma.NAD(b,:) = SetNormDiff( SomaBaseTestData, SomaRunTestData );
            ProcBaseTestData = periBout.proc.A(ZtestBase,:,b);
            ProcRunTestData = periBout.proc.A(ZtestRun,:,b);
            periBout.proc.NAD(b,:) = SetNormDiff( ProcBaseTestData, ProcRunTestData );
            % Bootstrap to estimate noise
            [~,baseBootSamp] =  bootstrp( Nboot, [], 1:NtestBase ); % resample baseline data
            [~,runBootSamp] =  bootstrp( Nboot, [], 1:NtestRun ); % resample run-onset data
            for n = flip(1:Nboot)
                CellBootDiff(n,:) = SetNormDiff( CellBaseTestData(baseBootSamp(:,n),:), CellRunTestData(runBootSamp(:,n),:) );   %#ok<*AGROW>
                SomaBootDiff(n,:) = SetNormDiff( SomaBaseTestData(baseBootSamp(:,n),:), SomaRunTestData(runBootSamp(:,n),:) );  
                ProcBootDiff(n,:) = SetNormDiff( ProcBaseTestData(baseBootSamp(:,n),:), ProcRunTestData(runBootSamp(:,n),:) );  
            end
            periBout.cell.NADci(b,:,2) = prctile( CellBootDiff, 90, 1 );  periBout.cell.NADci(b,:,1) = prctile( CellBootDiff, 10, 1 ); % Calculate confidence intervals
            periBout.cell.trig(b, periBout.cell.NADci(b,:,1) > 0) = true; % consider it a locomotion-triggered event if 90 percent CI excludes 0 difference
            periBout.cell.trigFrac(b) = sum(periBout.cell.trig(b,:), 2)/Ncell;
            periBout.soma.NADci(b,:,2) = prctile( SomaBootDiff, 90, 1 );  periBout.soma.NADci(b,:,1) = prctile( SomaBootDiff, 10, 1 ); 
            periBout.soma.trig(b, periBout.soma.NADci(b,:,1) > 0) = true; 
            periBout.soma.trigFrac(b) = sum(periBout.soma.trig(b,:), 2)/Ncell;
            periBout.proc.NADci(b,:,2) = prctile( ProcBootDiff, 90, 1 );  periBout.proc.NADci(b,:,1) = prctile( ProcBootDiff, 10, 1 ); 
            periBout.proc.trig(b, periBout.proc.NADci(b,:,1) > 0) = true; 
            periBout.proc.trigFrac(b) = sum(periBout.proc.trig(b,:), 2)/Ncell;
            % Quantify difference between peak activity during bout vs just before
            [cellResp, ~] = max(periBout.cell.A(goodTestRun,:,b), [], 1); % cellRespInd
            periBout.cell.delta(b,:) = cellResp - median( periBout.cell.A(goodTestBase,:,b), 1, 'omitnan' ); 
            periBout.cell.peak(b,:) = cellResp;
            [somaResp, ~] = max(periBout.soma.A(goodTestRun,:,b), [], 1);
            periBout.soma.peak(b,:) = somaResp;
            periBout.soma.delta(b,:) = somaResp - median( periBout.soma.A(goodTestBase,:,b), 1, 'omitnan' );
            [procResp, ~] = max(periBout.proc.A(goodTestRun,:,b), [], 1);
            periBout.proc.peak(b,:) = procResp;
            periBout.proc.delta(b,:) = procResp - median( periBout.proc.A(goodTestBase,:,b), 1, 'omitnan' );
            % Determine latency from locomotion onset to threshold crossing
            periBout.cell.thresh(b,:) = median( periBout.cell.A(goodTestBase,:,b)) + 2*mad( periBout.cell.A(goodTestBase,:,b), 1, 1); % triggering threshold = median base activity + 2 MAD
            periBout.soma.thresh(b,:) = median( periBout.soma.A(goodTestBase,:,b)) + 2*mad( periBout.soma.A(goodTestBase,:,b), 1, 1);
            periBout.proc.thresh(b,:) = median( periBout.proc.A(goodTestBase,:,b)) + 2*mad( periBout.proc.A(goodTestBase,:,b), 1, 1);
            for c = 1:Ncell 
                cellTrigFrame = goodTestRun(find( periBout.cell.A(goodTestRun,c,b) - periBout.cell.thresh(b,c) > 0, 1, 'first' ));
                if periBout.cell.trig(b,c) && ~isempty(cellTrigFrame), periBout.cell.latency(b,c) = periBout.T(cellTrigFrame ); end
                somaTrigFrame = goodTestRun(find( periBout.soma.A(goodTestRun,c,b) - periBout.soma.thresh(b,c) > 0, 1, 'first' ));
                if periBout.soma.trig(b,c) && ~isempty(somaTrigFrame), periBout.soma.latency(b,c) = periBout.T(somaTrigFrame ); end
                procTrigFrame = goodTestRun(find( periBout.proc.A(goodTestRun,c,b) - periBout.proc.thresh(b,c) > 0, 1, 'first' ));
                if periBout.proc.trig(b,c) && ~isempty(procTrigFrame), periBout.proc.latency(b,c) = periBout.T(procTrigFrame ); end
            end
        else
            fprintf('\nb = %i: not enough good test frames', b)
        end
        % Inspect the bout and fluor (optional)
        if inspectBout
            % Velocity
            subplot(sp(4)); 
            plot( periBout.T, periBout.velocity(:,b) ); hold on;      
            set(gca,'Xtick',-baseDur:1:runDur, 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
            Ylim = get(gca,'Ylim');
            line(testLim(1)*[1,1], Ylim, 'lineStyle','--', 'Color','r');
            line(testLim(2)*[1,1], Ylim, 'lineStyle','--', 'Color','r');
            hold off;
            ylabel('Velocity (cm/s)'); xlabel('Time Running (s)'); 
            title( sprintf('initMax = %2.2f cm/s', periBout.initMax(b) ) );
            for c = 1:Ncell
                % Soma activity
                subplot(sp(1));
                plot( periBout.T, periBout.soma.A(:,c,b) ); hold on;
                line([-10,10], periBout.soma.thresh(b,c)*[1,1], 'lineStyle','--', 'color','k'); 
                hold off;
                set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
                ylabel('Somatic Activity'); 
                title( sprintf('b = %i: soma %i. NAD = %2.1f (CI = [%2.1f, %2.1f]). trigger = %i. Latency = %2.1f s,  delta = %2.1f', ...
                    b, c, periBout.soma.NAD(b,c), periBout.soma.NADci(b,c,1), periBout.soma.NADci(b,c,2), periBout.soma.trig(b,c), periBout.soma.latency(b,c), periBout.soma.delta(b,c) )); %  
                % Process activity
                subplot(sp(2));
                plot( periBout.T, periBout.proc.A(:,c,b) ); hold on;
                line([-10,10], periBout.proc.thresh(b,c)*[1,1], 'lineStyle','--', 'color','k'); hold off;
                set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
                ylabel('Process Activity'); 
                title( sprintf('NAD = %2.1f (CI = [%2.1f, %2.1f]). trigger = %i. Latency = %2.1f s,  delta = %2.1f', ...
                     periBout.proc.NAD(b,c), periBout.proc.NADci(b,c,1), periBout.proc.NADci(b,c,2), periBout.proc.trig(b,c), periBout.proc.latency(b,c), periBout.proc.delta(b,c) )); %  
                % Cell activity
                subplot(sp(3)); 
                plot( periBout.T, periBout.cell.A(:,c,b) );  hold on;
                line([-10,10], periBout.cell.thresh(b,c)*[1,1], 'lineStyle','--', 'color','k'); hold off;
                set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
                ylabel('Cell Activity'); 
                title( sprintf('NAD = %2.1f (CI = [%2.1f, %2.1f]). trigger = %i. Latency = %2.1f s,  delta = %2.1f', ...
                     periBout.cell.NAD(b,c), periBout.cell.NADci(b,c,1), periBout.cell.NADci(b,c,2), periBout.cell.trig(b,c), periBout.cell.latency(b,c), periBout.cell.delta(b,c) )); %  
                linkaxes(sp,'x');
                xlim([min(periBout.T(1),-baseDur), max(periBout.T(end),runDur)]);
                pause;
            end           
        end
    end
    periInd.base = periBout.frames(1:periBout.Nbase,:); periInd.base = unique(periInd.base(~isnan(periInd.base)));
    periInd.run = periBout.frames(periBout.Nbase+1:end,:); periInd.run = unique(periInd.run(~isnan(periInd.run)));
    periInd.all = unique( periBout.frames(~isnan(periBout.frames)) ); %  unique( reshape( [periBout.frames], periBout.Nbout*periBout.Nframe, 1 ) );
    if show
        FS = 16;
        figure('WindowState','max', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        sp(1) = subplot(3,1,1);
        meanProc = mean(periBout.proc.A, 3, 'omitnan'); % average over all bouts
        plot( periBout.T, meanProc ); hold on;
        plot( periBout.T, mean(meanProc,2), 'k', 'LineWidth',2 ); % average over all bouts and ROI
        set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
        xlim([min(periBout.T(1), -baseDur), max(periBout.T(end), runDur)]); % xlim([-baseDur, runDur]); %
        ylabel('Mean Proc Fluor'); % xlabel('Time Running (s)'); 
        title( sprintf('%i cells', Ncell ) );

        sp(2) = subplot(3,1,2);
        meanSoma = mean(periBout.soma.A, 3, 'omitnan'); % average over all bouts
        plot( periBout.T, meanSoma ); hold on;
        plot( periBout.T, mean(meanSoma,2), 'k', 'LineWidth',2 ); % average over all bouts and ROI
        set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
        ylabel('Mean Soma Fluor'); % xlabel('Time Running (s)'); 
        xlim([min(periBout.T(1), -baseDur), max(periBout.T(end), runDur)]);

        sp(3) = subplot(3,1,3);
        plot( periBout.T, periBout.velocity, 'LineWidth',2 ); % 'Color',[0,0,0,AlphaVal] periSpeed
        xlim([min(periBout.T(1), -baseDur), max(periBout.T(end), runDur)]); 
        set(gca,'Xtick',-baseDur:2:runDur, 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
        ylabel('Velocity (cm/s)'); xlabel('Time Running (s)'); 
        title( sprintf('%i bouts', Nbout ) );
        linkaxes(sp,'x');
    end
end
end