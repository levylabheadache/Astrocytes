% Gather all the locomotive bouts and peribouts
baseDur = 10; runDur = 10;
locoBout = cell(1,Nexpt); Nbout = cell(1,Nexpt); periBout = cell(1,Nexpt);  periFrame = cell(1,Nexpt); 
for x = 1:Nexpt
    for f = flip(1:Nfov(x)) %fov{X}
        for m = 1:NfovMovies{x}(f) % mExpt{x}{f}
            [locoBout{x}{f}{m}, Nbout{x}{f}(m)] = LocoBouts(Tadj{x}{f}{m}, velocity{x}{f}{m} ); % , 'show',true, 'pause',true 
            [periBout{x}{f}(m), periFrame{x}{f}{m}] = GetPeriBout( Tadj{x}{f}{m}, velocity{x}{f}{m}, locoBout{x}{f}{m}, dFF{x}{f}(m), A{x}{f}(m), ...
                'base',baseDur, 'run',runDur, 'inspect',false, 'show',false); % 
        end
    end
    %[shortPeri{x}, shortBout{x}] = FilterPoolPeriBout(locoBout{x}, periBout{x}, mType{x}, shortPeriFilt, 'show',false ); 
end
shortDur = [0,1.5]; shortPeriFilt = struct('iso',[12,10], 'dur',shortDur);
[shortPeri, shortBout] = cellfun( @FilterPoolPeriBout, locoBout, periBout, mType, repmat({shortPeriFilt}, 1, Nexpt), 'UniformOutput',false ); 

%% Compare properties of brief running bouts to resulting (process) fluor events
varNames = {'MaxSpd', 'RunDur', 'Travel', 'TrigFrac', 'Delta', 'Lat'};
Nvar = numel(varNames);
shortPeriTable = array2table(zeros(0,Nvar), 'VariableNames',varNames);
r = 0;
show = false; % true; % 
if show, figure('WindowState','max'); end
for x = 1:Nexpt
    for f = find( [shortPeri{x}.Nbout] )
        % Exclude drug-affected 
        if contains( group{x}, 'Saline' ) || contains( group{x}, 'Ctrl' )
            bGet = [shortPeri{x}(f).type.pre, shortPeri{x}(f).type.stim];
        else
            bGet = shortPeri{x}(f).type.pre;
        end
        for b = bGet %[shortPeri{x}(f).type.pre, shortPeri{x}(f).type.stim]
            r = r + 1;
            % Bout properties
            shortPeriTable.MaxSpd(r) = shortPeri{x}(f).maxSpd(b); %shortBout{x}{f}(b).maxSpd; 
            shortPeriTable.RunDur(r) = shortPeri{x}(f).dur(b); % shortBout{x}{f}(b).dur; 
            shortPeriTable.Travel(r) = shortPeri{x}(f).travel(b); %shortBout{x}{f}(b).travel; 
            % Fluor response properties (averaged over all cells)
            shortPeriTable.TrigFrac(r) = shortPeri{x}(f).proc.trigFrac(b);
            shortPeriTable.Lat(r) = mean( shortPeri{x}(f).proc.latency(b,:), 2, 'omitnan' ); %  - shortBout{x}{f}(b).dur
            shortPeriTable.Delta(r) = mean( shortPeri{x}(f).proc.delta(b,:), 2, 'omitnan' ); %mean( shortPeri{x}(f).proc.delta(b,:), 2 );     
            if show
                subplot(2,1,2); 
                plot( shortPeri{x}(f).T, shortPeri{x}(f).velocity(:,b) ); xlim( shortPeri{x}(f).T([1,end]) );
                ylabel('Velocity (cm/s)'); xlabel('Peri-run Time (s)');
                title( sprintf('MaxSpd = %2.2f cm/s, RunDur = %2.2f s, Travel = %2.2f cm', shortPeriTable.MaxSpd(r), shortPeriTable.RunDur(r), shortPeriTable.Travel(r)  ) );
                set(gca,'Xtick',-20:2:20);
                subplot(2,1,1);
                if sum(shortPeri{x}(f).trig(b,:)), plot( shortPeri{x}(f).T, shortPeri{x}(f).proc.F(:,shortPeri{x}(f).trig(b,:),b), 'b' ); hold on; end
                if sum(~shortPeri{x}(f).trig(b,:)), plot( shortPeri{x}(f).T, shortPeri{x}(f).proc.F(:,~shortPeri{x}(f).trig(b,:),b), 'r' ); hold on; end
                title( sprintf('r = %i:  TrigFrac = %1.2f, TrigLat = %2.2f s, TrigMax = %2.2f', r, shortPeriTable.TrigFrac(r), shortPeriTable.TrigLat(r), shortPeriTable.TrigMax(r)  ) );
                h(1) = plot(1,NaN, 'b'); h(2) = plot(1,NaN, 'r'); 
                legend(h, {'Run-triggered', 'Non-triggered'}, 'Location','NorthWest');
                set(gca,'Xtick',-20:2:20);
                hold off;
                xlim( shortPeri{x}(f).T([1,end]) );
                pause;
            end
        end
    end
end
Nshort = size(shortPeriTable,1);
fprintf('\n\nFound %i well-isolated short-duration bouts', Nshort);

%% Plot velocity, and somatic and process fluor for each short bout
clearvars h sp;
figure('WindowState','max'); 
sp(3) = subplot(3,1,3); sp(2) = subplot(3,1,2); sp(1) = subplot(3,1,1);
for x = 1:Nexpt
    for f = find( [shortPeri{x}.Nbout] )
        for b = [shortPeri{x}(f).type.pre, shortPeri{x}(f).type.stim]
            subplot(sp(1));
            Ntrig = sum(  shortPeri{x}(f).proc.trig(b,:) );
            if Ntrig == 0
                plot( shortPeri{x}(f).T, shortPeri{x}(f).proc.A(:,~shortPeri{x}(f).proc.trig(b,:),b), 'LineWidth',1, 'color','r'); hold on;
                plot( shortPeri{x}(f).T, mean(shortPeri{x}(f).proc.A(:,~shortPeri{x}(f).proc.trig(b,:),b),2), 'LineWidth',2, 'Color','r' ); 
            elseif Ntrig == Ncell{x}(f)
                plot( shortPeri{x}(f).T, shortPeri{x}(f).proc.A(:,shortPeri{x}(f).proc.trig(b,:),b), 'LineWidth',1, 'color','b' ); hold on;
                plot( shortPeri{x}(f).T, mean(shortPeri{x}(f).proc.A(:,shortPeri{x}(f).proc.trig(b,:),b),2), 'LineWidth',2, 'Color','b' ); 
            else
                plot( shortPeri{x}(f).T, shortPeri{x}(f).proc.A(:,~shortPeri{x}(f).proc.trig(b,:),b), 'LineWidth',1, 'color','r' ); hold on;
                plot( shortPeri{x}(f).T, shortPeri{x}(f).proc.A(:,shortPeri{x}(f).proc.trig(b,:),b), 'LineWidth',1, 'color','b' ); 
                plot( shortPeri{x}(f).T, mean(shortPeri{x}(f).proc.A(:,~shortPeri{x}(f).proc.trig(b,:),b),2), 'LineWidth',2, 'Color','r' ); 
                plot( shortPeri{x}(f).T, mean(shortPeri{x}(f).proc.A(:,shortPeri{x}(f).proc.trig(b,:),b),2), 'LineWidth',2, 'Color','b' ); 
            end

            ylabel('Fproc'); xlim(shortPeri{x}(f).T([1,end]));
            hold off;

            subplot(sp(2)); 
            plot( shortPeri{x}(f).T, shortPeri{x}(f).soma.A(:,:,b), 'LineWidth',1 ); hold on;
            plot( shortPeri{x}(f).T, mean(shortPeri{x}(f).soma.A(:,:,b),2), 'LineWidth',2, 'Color','k' );             
            ylabel('Fsoma'); xlim(shortPeri{x}(f).T([1,end]));
            hold off;
            
            subplot(sp(3));
            plot( shortPeri{x}(f).T, shortPeri{x}(f).velocity(:,b), 'LineWidth',1.5, 'Color','k' ); %hold on;  % [0,0,0,setAlpha]
            xlabel('Peri-bout time'); ylabel('Velocity (cm/s)'); xlim(shortPeri{x}(f).T([1,end]));
            pause;
        end
    end
end


%%  Plot all filtered periBouts (velocity and fluor) - pre-and post-stimulus separated
clearvars h sp;
figure('WindowState','max'); 
SP(2) = subplot(2,2,2); SP(1) = subplot(2,2,1); sp(2) = subplot(2,2,4); sp(1) = subplot(2,2,3);
setAlpha = 0.2;
for x = 1:Nexpt
    for f = find( [shortPeri{x}.Nbout] )
        for b = shortPeri{x}(f).type.pre
            subplot(sp(1)); 
            plot( shortPeri{x}(f).T, shortPeri{x}(f).velocity(:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
            subplot(SP(1));
            plot( shortPeri{x}(f).T, shortPeri{x}(f).cell.A(:,:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
        end
        subplot(sp(1)); ylabel('Velocity (cm/s)'); xlim([-baseDur, runDur]);
        subplot(sp(2)); ylabel('Velocity (cm/s)'); xlim([-baseDur, runDur]);
        
        for b = shortPeri{x}(f).type.stim
            subplot(sp(2));
            plot( shortPeri{x}(f).T, shortPeri{x}(f).velocity(:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
            subplot(SP(2));
            plot( shortPeri{x}(f).T, shortPeri{x}(f).cell.A(:,:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
        end
        subplot(SP(1)); ylabel('Activity'); xlim([-baseDur, runDur]); title( sprintf('%s %s FOV %i: %s', mouse{x}, exptDate{x}, fov{x}(f), group{x} ) ); 
        subplot(SP(2)); ylabel('Activity'); xlim([-baseDur, runDur]);
        
        pause;
        subplot(sp(1)); cla;  subplot(sp(2)); cla;  subplot(SP(1)); cla; subplot(SP(2)); cla; 
    end
end

%% Which subROI (soma or procees) is most sensitive to locomotion?
shortPeri{x}(f).soma.trigFrac



%%
clearvars h sp;
figure('WindowState','max'); 
% Multivariate linear regression -> TrigFrac 
%mvFit = fitlm( [shortPeriTable.MaxSpd, shortPeriTable.RunDur, shortPeriTable.Travel ], shortPeriTable.TrigFrac );

% MaxSpd vs TrigFrac
subplot(3,3,1);
plot( shortPeriTable.MaxSpd, shortPeriTable.TrigFrac, '.' ); hold on;
axis square; box off;
xlabel('Peak Speed (cm/s)'); ylabel('Fraction of Cells Triggered');
ylim([0,1.01]);
% single-variable linear regression
% {
tempFit = fitlm( shortPeriTable.MaxSpd, shortPeriTable.TrigFrac );
xMax = max(shortPeriTable.MaxSpd); xMin = min(shortPeriTable.MaxSpd);
line( [xMin,xMax], predict(tempFit, [xMin;xMax]), 'Color','k', 'LineStyle','--'  ) 
title( sprintf('Short Bouts (duration %1.0f - %1.0f s, n=%i). R^2 = %2.2f', shortDur(1), shortDur(2), Nshort, tempFit.Rsquared.Ordinary));
%}
% RunDur vs TrigFrac
subplot(3,3,2);
plot( shortPeriTable.RunDur, shortPeriTable.TrigFrac, '.' );
axis square; box off;
xlabel('Run Duration (s)'); ylabel('Fraction of Cells Triggered');
ylim([0,1.01]); %xlim([0,Inf]);
tempFit = fitlm( shortPeriTable.RunDur, shortPeriTable.TrigFrac );
xMax = max(shortPeriTable.RunDur); xMin = min(shortPeriTable.RunDur);
line( [xMin,xMax], predict(tempFit, [xMin;xMax]), 'Color','k', 'LineStyle','--'  ) 
title( sprintf('R^2 = %2.2f', tempFit.Rsquared.Ordinary));
 
% Travel vs TrigFrac
subplot(3,3,3);
plot( shortPeriTable.Travel, shortPeriTable.TrigFrac, '.' );
axis square; box off;
xlabel('Travel Distance (cm)'); ylabel('Fraction of Cells Triggered');
ylim([0,1.01]);
tempFit = fitlm( shortPeriTable.Travel, shortPeriTable.TrigFrac );
xMax = max(shortPeriTable.Travel); xMin = min(shortPeriTable.Travel);
line( [xMin,xMax], predict(tempFit, [xMin;xMax]), 'Color','k', 'LineStyle','--'  ) 
title( sprintf('R^2 = %2.2f', tempFit.Rsquared.Ordinary));

% MaxSpd vs mean delta dF/F (peak vs baseline)
subplot(3,3,4);
plot( shortPeriTable.MaxSpd, shortPeriTable.Delta, '.' );
axis square; box off;
xlabel('Peak Speed (cm/s)'); ylabel('Delta dF/F');
% RunDur vs mean delta dF/F (peak vs baseline)
subplot(3,3,5);
plot( shortPeriTable.RunDur, shortPeriTable.Delta, '.' );
axis square; box off;
xlabel('Run Duration (s)'); ylabel('Delta dF/F');
% Travel vs mean delta dF/F (peak vs baseline)
subplot(3,3,6);
plot( shortPeriTable.Travel, shortPeriTable.Delta, '.' );
axis square; box off;
xlabel('Travel Distance (cm)'); ylabel('Delta dF/F');

% MaxSpd vs Latency
subplot(3,3,7);
plot( shortPeriTable.MaxSpd, shortPeriTable.Lat, '.' );
axis square; box off;
xlabel('Peak Speed (cm/s)'); ylabel('Latency to Trigger (s)');

% RunDur vs Latency
subplot(3,3,8);
plot( shortPeriTable.RunDur, shortPeriTable.Lat, '.' );
axis square; box off;
xlabel('Run Duration (s)'); ylabel('Latency to Trigger (s)');

% Travel vs Latency
subplot(3,3,9);
plot( shortPeriTable.Travel, shortPeriTable.Lat, '.' );
axis square; box off;
xlabel('Travel Distance (cm)'); ylabel('Latency to Trigger (s)');