%% Compare locomotion/peribout responses before/after injection
baseBouts = cell(1,Nexpt); stimBouts = cell(1,Nexpt); %[];
baseTravelTot = nan(1,Nexpt); stimTravelTot = nan(1,Nexpt);
baseCellNADmean = cell(1,Nexpt); stimCellNADmean = cell(1,Nexpt);
for g = 1:Ngroup
    for x = xGroup{g}
        for f = flip(1:Nfov(x)) %
            %
            baseBouts{x} = [baseBouts{x}, [locoBout{x}{f}{mType{x}(f).base}]];
            stimBouts{x} = [stimBouts{x}, [locoBout{x}{f}{mType{x}(f).stim}]];
            % Gather isolated periBouts
            baseCellNADmean{x}{f} = mean( isoPeri{x}(f).cell.NAD(isoPeri{x}(f).type.pre,:), 1, 'omitnan' );
            stimCellNADmean{x}{f} = mean( isoPeri{x}(f).cell.NAD(isoPeri{x}(f).type.stim,:), 1, 'omitnan' );
        end
        if ~isempty(baseBouts{x})
            baseTravelTot(x) = sum([baseBouts{x}.travel]); 
        elseif baseDurTot(x) > 0
            baseTravelTot(x) = 0;
        end
        if ~isempty(stimBouts{x}) 
            stimTravelTot(x) = sum([stimBouts{x}.travel]); 
        elseif stimDurTot(x) > 0
            stimTravelTot(x) = 0;
        end
        % Compare properties of locoBouts
    end
end
baseSpeedAvg = baseTravelTot./baseDurTot;
stimSpeedAvg = stimTravelTot./stimDurTot;
%% Average speed before/after stimulus
figure('WindowState','maximized');
for g = flip(1:Ngroup)
    subplot(1,Ngroup, g);
    plot( [1,2], [baseSpeedAvg(xGroup{g})', stimSpeedAvg(xGroup{g})'] ); hold on;
    plot( [1,2], [baseSpeedAvg(xGroup{g})', stimSpeedAvg(xGroup{g})'], '.' );
    title( sprintf('%s', setGroup{g} ) );
    axis square;
    xlim([0.5,2.5]); ylim([0,10]);
    set(gca,'Xtick',1:2, 'XtickLabel',{'Base','Stim'}, 'TickDir','out');
end
ylabel( 'Mean Speed (cm/s)');

%% Average NAD before/after stimulus
figure('WindowState','maximized');
for g = flip(1:Ngroup)
    for x = xGroup{g}
        for f = 1:Nfov(x) %find( [isoPeri{x}.Nbout] )
            subplot(1,Ngroup, g);
            plot( [1,2], [baseCellNADmean{x}{f}', stimCellNADmean{x}{f}'] ); hold on;
            plot( [1,2], [baseCellNADmean{x}{f}', stimCellNADmean{x}{f}'], '.' );
            title( sprintf('%s', setGroup{g} ) );
            axis square;
            xlim([0.5,2.5]); ylim([-1,1]);
            set(gca,'Xtick',1:2, 'XtickLabel',{'Base','Stim'}, 'TickDir','out');
        end
    end
end
ylabel( 'Mean Normalized Activity Difference (cm/s)');

%% Individual periBouts before/after stimulus
for g = 2:Ngroup
    for x = xGroup{g}
        for f = flip(1:Nfov(x)) %
            for r = 1:Ncell{x}(f)
                % Plot isolated periBouts
                if ~isempty( isoPeri{x}(f).type.pre )
                    plot( isoPeri{x}(f).T, squeeze(isoPeri{x}(f).cell.A(:,r,isoPeri{x}(f).type.pre)), 'k' ); hold on;
                end
                if ~isempty( isoPeri{x}(f).type.stim )
                    plot( isoPeri{x}(f).T, squeeze(isoPeri{x}(f).cell.A(:,r,isoPeri{x}(f).type.stim)), 'r' );
                end
                title( sprintf('%s %s FoV %i cell %i', mouse{x}, exptDate{x}, fov{x}(f), r ) );
                xlim([-10, 10]);
                pause;
                cla;
            end
        end
    end
end