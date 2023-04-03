function [proj, cellROI, somaROI, procROI, backROI] = ReviewAstroROI(mouse, date, fov, varargin) %#ok<*STOUT>


checkFOV = @(x)(ischar(x) || isnumeric(x));
IP = inputParser;
addRequired( IP, 'mouse', @ischar) 
addRequired( IP, 'date', @ischar) 
addRequired( IP, 'fov', checkFOV) 
addParameter( IP, 'show', true, @islogical )
addParameter( IP, 'add', false, @islogical )
parse( IP, mouse, date, fov, varargin{:} );  
addROI = IP.Results.add; 
show = IP.Results.show; 
if addROI, show = true; end

[dataName, dataDir] = SetDataPath( mouse, date, [], fov);
roiMat = FileFind( dataDir, 'mat', false, @(x)(contains( x, 'ROI' )) );
if ~isempty( roiMat )
    load(roiMat{2} );
    if show
        opt = {[0.06,0.03], [0.04,0.02], [0.01,0.01]};  % {[vert, horz], [bottom, top], [left, right] } 
        FS = 12;
        ColorMat = distinguishable_colors(Ncell);
        close all; clearvars sp h;
        figure('units','normalized','WindowState','max', 'color','k'); % , 'name',movieName
        sp(1) = subtightplot(1,3,1,opt{:});
        imshow( proj.mean, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
        plot( backROI.edge(:,1), backROI.edge(:,2), '--', 'Color', 'w' );  
        title('Mean Projection', 'FontSize',FS, 'Color','w');
        for c = 1:Ncell
            for n = 1:somaROI(c).Nedge %#ok<*NODEF>
                plot( somaROI(c).edge{n}(:,1), somaROI(c).edge{n}(:,2), '--', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end  
            for n = 1:cellROI(c).Nedge
                plot( cellROI(c).edge{n}(:,1), cellROI(c).edge{n}(:,2), '-.', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end
        end

        sp(2) = subtightplot(1,3,2,opt{:});
        imshow( proj.max, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
        plot( backROI.edge(:,1), backROI.edge(:,2), '--', 'Color', 'w' );  
        for c = 1:Ncell
            for n = 1:somaROI(c).Nedge
                plot( somaROI(c).edge{n}(:,1), somaROI(c).edge{n}(:,2), '--', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end  
            for n = 1:cellROI(c).Nedge
                plot( cellROI(c).edge{n}(:,1), cellROI(c).edge{n}(:,2), '-.', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end
        end

        title('Max Projection', 'FontSize',FS, 'Color','w');
        sp(3) = subtightplot(1,3,3,opt{:});
        imshow( proj.std, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
        plot( backROI.edge(:,1), backROI.edge(:,2), '--', 'Color', 'w' );  
        for c = 1:Ncell
            for n = 1:somaROI(c).Nedge
                plot( somaROI(c).edge{n}(:,1), somaROI(c).edge{n}(:,2), '--', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end  
            for n = 1:cellROI(c).Nedge
                plot( cellROI(c).edge{n}(:,1), cellROI(c).edge{n}(:,2), '-.', 'Color', ColorMat(c,:) );  % colorblind(c,:)
            end
        end
        title('Standard Deviation Projection', 'FontSize',FS, 'Color','w');
        linkaxes(sp,'xy');
        if addROI
            [newCell, newSoma, newProc, Nnew] = DrawAstroROI();
            Ncell = Ncell + Nnew;
            cellROI = [cellROI, newCell];
            somaROI = [somaROI, newSoma];
            procROI = [procROI, newProc];
            % Write cells to temporary variables for saving
            savePath = sprintf('%s%s_ROI.mat', dataDir, dataName); % metadata(1).root
            if exist('roiParam','var')
                save(savePath, 'backROI', 'cellROI', 'somaROI', 'procROI', 'Ncell', 'proj', 'roiParam' ); 
            else
                save(savePath, 'backROI', 'cellROI', 'somaROI', 'procROI', 'Ncell', 'proj' ); 
            end
            fprintf('\n  Saved %s\n',savePath); 
            roi2aqua( mouse, date, fov );
        end
    end
else
    error('No ROI file to load!');
end
end