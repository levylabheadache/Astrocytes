% Gather all bouts
allBouts = [];
for x = 1:Nexpt
    for f = 1:Nfov(x)
        for m = find(Nbout{x}{f} > 0)
            allBouts = [allBouts, locoBout{x}{f}{m}];
        end
    end
end
NallBouts = numel(allBouts);
%% Distributions of properties of locomotive bouts
durBin = 1;
travBin = 10;
spdBin = 1;
isoBin = 2;

figure('WindowState', 'maximized');
subplot(2,2,1);
histogram( [allBouts.dur], 'BinWidth', durBin, 'Normalization', 'probability'  );
set(gca, 'Xtick', 0:5:50, 'TickDir','out', 'box','off' ); axis square;
xlabel('Duration (s)'); title( sprintf('n = %i bouts', NallBouts) );
xlim([0,40]);
subplot(2,2,2);
histogram( [allBouts.travel], 'BinWidth', travBin, 'Normalization', 'probability'  );
set(gca, 'Xtick', 0:50:500, 'TickDir','out', 'box','off' );
xlabel('Distance Traveled (cm)'); axis square;
xlim([0,300]);
subplot(2,2,3);
histogram( [allBouts.maxSpd], 'BinWidth', spdBin, 'Normalization', 'probability' );
set(gca, 'Xtick', 0:5:50, 'TickDir','out', 'box','off' );
xlabel('Top Speed (cm/s)'); axis square;
xlim([0,50]);
subplot(2,2,4);
allIso = reshape([allBouts.iso], NallBouts, 2);
histogram( allIso(:,1), 'BinWidth', isoBin, 'Normalization', 'cdf'  );
set(gca, 'Xtick', 0:10:120, 'TickDir','out', 'box','off' );
xlabel('Isolation (s)'); axis square;
xlim([0,120]);

%% Compare well-isolated bouts that trigger fluor events to non-triggering
isoPeriFilt = struct('iso',[10,10], 'dur',[0,Inf]); % 
[isoPeri, ~] = cellfun( @FilterPoolPeriBout, locoBout, periBout, mType, repmat({isoPeriFilt}, 1, Nexpt), 'UniformOutput',false );  
% Pool properties of brief running bouts to resulting (process) fluor events
varNames = {'MaxSpd', 'RunDur', 'Travel', 'SomaFrac', 'ProcFrac'}; % , 'TrigFrac', 'Delta', 'Lat'
Nvar = numel(varNames);
isoPeriTable = array2table(zeros(0,Nvar), 'VariableNames',varNames);
r = 0;
for x = 1:Nexpt
    for f = find( [isoPeri{x}.Nbout] )
        % Exclude drug-affected 
        if contains( group{x}, 'Saline' ) || contains( group{x}, 'Ctrl' )
            bGet = [isoPeri{x}(f).type.pre, isoPeri{x}(f).type.stim];
        else
            bGet = isoPeri{x}(f).type.pre;
        end
        for b = bGet %[isoPeri{x}(f).type.pre, isoPeri{x}(f).type.stim]
            r = r + 1;
            % Bout properties
            isoPeriTable.MaxSpd(r) = isoPeri{x}(f).maxSpd(b); 
            isoPeriTable.RunDur(r) = isoPeri{x}(f).dur(b);  
            isoPeriTable.Travel(r) = isoPeri{x}(f).travel(b); 
            % Fluor response properties (averaged over all cells)
            isoPeriTable.SomaFrac(r) = isoPeri{x}(f).soma.trigFrac(b);
            isoPeriTable.ProcFrac(r) = isoPeri{x}(f).proc.trigFrac(b); 
        end
    end
end
Niso = size(isoPeriTable,1);
fprintf('\n\nFound %i well-isolated bouts', Niso);

%% How do bout properties correlate with triggering (soma or process)
close all; clearvars h sp;
figure('WindowState','max'); 
% Soma
subplot(2,3,1);
plot( isoPeriTable.RunDur, isoPeriTable.SomaFrac, '.' );
xlabel('Bout Duration (s)');  ylabel('Somatic Triggering Fraction');
subplot(2,3,2);
plot( isoPeriTable.MaxSpd, isoPeriTable.SomaFrac, '.' );
xlabel('Max Speed (cm/s)'); 
subplot(2,3,3);
plot( isoPeriTable.Travel, isoPeriTable.SomaFrac, '.' );
xlabel('Distance Traveled (cm)'); 
% Process
subplot(2,3,4);
plot( isoPeriTable.RunDur, isoPeriTable.ProcFrac, '.' );
xlabel('Bout Duration (s)');  ylabel('Process Triggering Fraction');
subplot(2,3,5);
plot( isoPeriTable.MaxSpd, isoPeriTable.ProcFrac, '.' );
xlabel('Max Speed (cm/s)'); 
subplot(2,3,6);
plot( isoPeriTable.Travel, isoPeriTable.ProcFrac, '.' );
xlabel('Distance Traveled (cm)');

%% Sigmoidal fits
close all
%StimResponse( isoPeriTable.RunDur, isoPeriTable.SomaFrac, 0, 10, 'show',true, 'fit',true, 'names',{'Duration (s)', 'SomaFrac'} );
StimResponse( isoPeriTable.RunDur, isoPeriTable.ProcFrac, 0, 10, 'show',true, 'fit',true, 'names',{'Duration (s)', 'ProcFrac'} );
%StimResponse( isoPeriTable.Travel, isoPeriTable.SomaFrac, 0, 10, 'show',true, 'fit',true, 'names',{'Travel (cm)', 'SomaFrac'} );
StimResponse( isoPeriTable.Travel, isoPeriTable.ProcFrac, 0, 10, 'show',true, 'fit',true, 'names',{'Travel (cm)', 'ProcFrac'} );
%StimResponse( isoPeriTable.MaxSpd, isoPeriTable.SomaFrac, 0, 10, 'show',true, 'fit',true, 'names',{'Max Speed (cm/s)', 'SomaFrac'} );
StimResponse( isoPeriTable.MaxSpd, isoPeriTable.ProcFrac, 0, 10, 'show',true, 'fit',true, 'names',{'Max Speed (cm/s)', 'ProcFrac'} );

%%
figure('WindowState','max'); 
for x = flip(1:Nexpt)
    for f = 1:Nfov(x)
        for b = 1:isoPeri{x}(f).Nbout
            subplot(2,1,2); 
            plot( isoPeri{x}(f).T, isoPeri{x}(f).velocity(:,b) ); xlim( isoPeri{x}(f).T([1,end]) );
            ylabel('Velocity (cm/s)'); xlabel('Peri-run Time (s)');
            set(gca,'Xtick',[-10:2:20]);
            subplot(2,1,1);
            if sum(isoPeri{x}(f).proc.trig(b,:)), plot( isoPeri{x}(f).T, isoPeri{x}(f).proc.A(:,isoPeri{x}(f).trig(b,:),b), 'b' ); hold on; end
            if sum(~isoPeri{x}(f).proc.trig(b,:)), plot( isoPeri{x}(f).T, isoPeri{x}(f).proc.A(:,~isoPeri{x}(f).trig(b,:),b), 'r' ); hold on; end
            h(1) = plot(1,NaN, 'b'); h(2) = plot(1,NaN, 'r'); 
            legend(h, {'Run-triggered', 'Non-triggered'}, 'Location','NorthWest');
            set(gca,'Xtick',[-10:2:20]);
            hold off;
            xlim( isoPeri{x}(f).T([1,end]) );
            pause;
        end
    end
end

%%  Plot all filtered periBouts (velocity and fluor) - pre-and post-stimulus separated
clearvars h sp;
figure('WindowState','max'); 
SP(2) = subplot(2,2,2); SP(1) = subplot(2,2,1); sp(2) = subplot(2,2,4); sp(1) = subplot(2,2,3);
setAlpha = 0.1;
for x = 1:Nexpt
    for f = find( [filtPeri{x}.Nbout] )
        for b = filtPeri{x}(f).type.pre
            subplot(sp(1)); 
            plot( filtPeri{x}(f).T, filtPeri{x}(f).velocity(:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
            subplot(SP(1));
            plot( filtPeri{x}(f).T, filtPeri{x}(f).cell(:,:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
        end
        for b = filtPeri{x}(f).type.stim
            subplot(sp(2));
            plot( filtPeri{x}(f).T, filtPeri{x}(f).velocity(:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
            subplot(SP(2));
            plot( filtPeri{x}(f).T, filtPeri{x}(f).cell(:,:,b), 'LineWidth',1.5, 'Color',[0,0,0,setAlpha] ); hold on; 
        end
        pause;
        subplot(sp(1)); cla;  subplot(sp(2)); cla;  subplot(SP(1)); cla; subplot(SP(2)); cla; 
    end
end
