function mergedPeriBouts = MergePeriBouts( locoBout, periBout, varargin )

IP = inputParser;
addRequired( IP, 'locoBout', @iscell )
addRequired( IP, 'periBout', @isstruct )
addParameter( IP, 'iso', 20, @isnumeric )
addParameter( IP, 'dur', [0,Inf], @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, locoBout, periBout, varargin{:} ); 
minIso = IP.Results.iso;
durRange = IP.Results.dur;
show = IP.Results.show; 

periBout(cellfun(@isempty, locoBout)) = [];
locoBout(cellfun(@isempty, locoBout)) = [];

totBout = sum([periBout.Nbout]);
Nperi = numel(periBout);
if Nperi > 0 % Nloco > 0 && 
    Ncell = size(periBout(1).cell,2);
    mergedPeriBouts = periBout(1);
    mergedPeriBouts.frames = nan(mergedPeriBouts.Nframe,0); 
    mergedPeriBouts.velocity = nan(mergedPeriBouts.Nframe,0); 
    mergedPeriBouts.cell = nan(mergedPeriBouts.Nframe,Ncell,0);  % nan(mergedPeriBouts.Nframe,Ncell,0); 
    mergedPeriBouts.soma = nan(mergedPeriBouts.Nframe,Ncell,0); 
    mergedPeriBouts.proc = nan(mergedPeriBouts.Nframe,Ncell,0);
    Ngood = zeros(1,Nperi);
    for p = 1:Nperi
        isoMat = vertcat(locoBout{p}.iso);
        isoMin = min(isoMat,[],2)';
        durVec = [locoBout{p}.dur];
        bGood = find( isoMin >= minIso & durVec >= durRange(1) & durVec <= durRange(2) );
        Ngood(p) = numel(bGood);
        mergedPeriBouts.velocity = cat(2, mergedPeriBouts.velocity, periBout(p).velocity(:,bGood) );
        mergedPeriBouts.cell = cat(3, mergedPeriBouts.cell, periBout(p).cell(:,:,bGood) );
        mergedPeriBouts.soma = cat(3, mergedPeriBouts.soma, periBout(p).soma(:,:,bGood) );
        mergedPeriBouts.proc = cat(3, mergedPeriBouts.proc, periBout(p).proc(:,:,bGood) ); 
    end
    totGood = sum(Ngood);
    fprintf('\nmerged %i of %i peri-locomotive bouts', totGood, totBout);

    if show && totGood > 0
        %close all;
        figure('WindowState','max', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        FS = 16;
        sp(1) = subplot(3,1,1);
        meanProc = mean(mergedPeriBouts.proc, 3, 'omitnan'); % average over all bouts
        plot( mergedPeriBouts.T, meanProc ); hold on;
        plot( mergedPeriBouts.T, mean(meanProc,2), 'k', 'LineWidth',2 ); % average over all bouts and ROI
        set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
        ylabel('Bout-Averaged Cell dF/F');
        title( sprintf('%i cells', Ncell ) );

        sp(2) = subplot(3,1,2);
        meanSoma = mean(mergedPeriBouts.soma, 3, 'omitnan'); % average over all bouts
        plot( mergedPeriBouts.T, meanSoma ); hold on;
        plot( mergedPeriBouts.T, mean(meanSoma,2), 'k', 'LineWidth',2 ); % average over all bouts and ROI
        set(gca,'Xtick',[], 'FontSize',FS, 'TickDir','out', 'box','off'); 
        ylabel('Bout-Averaged Somatic dF/F');

        sp(3) = subplot(3,1,3);
        plot( mergedPeriBouts.T, mergedPeriBouts.velocity ); 
        set(gca,'Xtick',round(mergedPeriBouts.T(1)):2:round(mergedPeriBouts.T(end)), 'FontSize',FS, 'TickDir','out', 'box','off'); % axis square;
        ylabel('Velocity (cm/s)'); xlabel('Time Running (s)'); 
        title( sprintf('%i bouts (%2.1f seconds isolation, duration %2.1f - %2.1f seconds)', totGood, minIso, durRange(1), durRange(2) ) );
        linkaxes(sp,'x');
        xlim(round( mergedPeriBouts.T([1,end]) )); 
    end
else
    mergedPeriBouts = struct('Nbout',NaN, 'Nframe',NaN, 'base',[], 'Nbase',NaN, 'run',[], 'Nrun',NaN, 'T',[], 'frames',[], 'velocity',[], 'cell',[], 'soma',[], 'proc',[]);
end

end

