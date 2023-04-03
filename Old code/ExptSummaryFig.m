%% time course: full fluor data plus locomotion
Tshow = cell(1,Nexpt); sp = cell(1,Nexpt); projSP = cell(1,Nexpt);
close all; clearvars sp projSP;
imOpt = {[0.01,0.01], [0.01,0.01], [0.01,0.01]};  % {[vert, horz], [bottom, top], [left, right] } 
traceOpt = {[0.06,0.03], [0.05,0.04], [0.04,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
FS = 12;
for x = flip(1:Nexpt)
    % Set up the figure
    figure('units','normalized','WindowState','max', 'color','w');
    for f = 1:Nfov(x)+1
        sp{x}(f) = subtightplot(Nfov(x)+1, 2, 2*f, traceOpt{:});
        projSP{x}(f) = subtightplot(Nfov(x)+1, 2, 2*f-1, imOpt{:});
    end
    % Plot the data
    for f = flip(1:Nfov(x))
        ColorMat = distinguishable_colors(Ncell{x}(f));
        [proj{x}(f), cellROI{x}{f}, somaROI{x}{f}, procROI{x}{f}, backROI{x}(f)] = ReviewAstroROI( mouse{x}, exptDate{x}, fov{x}(f), 'show',false, 'add',false );
        subplot(projSP{x}(f));
        imshow( proj{x}(f).max, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
        plot( backROI{x}(f).edge(:,1), backROI{x}(f).edge(:,2), '--', 'Color', 'w' );  
        for c = 1:Ncell{x}(f)
            for n = 1:somaROI{x}{f}(c).Nedge
                plot( somaROI{x}{f}(c).edge{n}(:,1), somaROI{x}{f}(c).edge{n}(:,2), '--', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end  
            for n = 1:cellROI{x}{f}(c).Nedge
                plot( cellROI{x}{f}(c).edge{n}(:,1), cellROI{x}{f}(c).edge{n}(:,2), '-.', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end
        end
        % Fluor
        Tshow{x}{f} = vertcat(Tmin{x}{f}{mType{x}(f).expt});
        subplot(sp{x}(f)) %subtightplot(Nfov(x)+1, 1, f, opt{:});
        plot( Tshow{x}{f}, vertcat(Fback{x}{f}{mType{x}(f).expt}), 'color',0.5*[1,1,1] ); hold on;
        plot( Tshow{x}{f}, vertcat(Fastro{x}{f}(mType{x}(f).expt).cell) );
        ylabel('Raw Fluorescence'); title( sprintf('FOV %i', fov{x}(f) ) );
        % Locomotion
        subplot(sp{x}(Nfov(x)+1))
        if ~any(cellfun( @isempty, velocity{x}{f} ))
            for m = mType{x}(f).expt, plot( Tmin{x}{f}{m}, velocity{x}{f}{m}, 'k' ); hold on;  end
        end
        %pause;
    end
    xlabel( 'Time After Injection (min)' ); ylabel( 'Velocity (cm/s)' ); 
    %linkaxes(projSP{x}, 'xy'); % axis tight;
    linkaxes(sp{x}, 'x'); axis tight;
    subplot(sp{x}(1));  
    xlim( [min(cellfun(@min, Tshow{x} )), max(cellfun(@max, Tshow{x} ))]);
    title( sprintf('%s %s FOV %i', mouse{x}, exptDate{x}, fov{x}(f)) ); % sprintf('%s %s: %2.2f mg/kg %s.  FOV %i', mouse{x}, exptDate{x}, dose{x}, drug{x}, fov{x}(f) )
    pause;
end
