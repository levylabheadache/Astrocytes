function [backROI, cellROI, somaROI, procROI, Ncell, proj] = MakeAstroROI( expt, varargin )
IP = inputParser;
addRequired( IP, 'expt', @isstruct)
addParameter( IP, 'conn', 8, @isnumeric )
addParameter( IP, 'radius', 5, @isnumeric ) %minRad = 5; % um
addParameter( IP, 'Nsigma', 3, @isnumeric ) % Nsigma = 2.5;
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, varargin{:} );
roiParam.conn = IP.Results.conn;
roiParam.radius = IP.Results.radius;
roiParam.Nsigma = IP.Results.Nsigma;
overwrite = IP.Results.overwrite;

% Convert pixels to microns
minPix = roiParam.radius/expt.umPerPixel; % express minimum radius (which is defined in um) in pixels
minArea = round( pi*(minPix)^2 );
roiPath = sprintf('%s%s_ROI.mat', expt.dir, expt.name);
projPath = sprintf('%sProjections\\%s_reg_green.tif', expt.dir, expt.name); %sprintf('%sProjections\\%s_greenProj_cat.tif', expt.dir, expt.name);
if exist(roiPath, 'file') && ~overwrite
    fprintf('\nLoading %s', roiPath)
    load(roiPath)
else
    % Generate projections and histograms
    greenMovie = loadtiff(projPath);
    movieRes = size(greenMovie, [1,2]);
    proj.mean = mean( greenMovie, 3 ); %metadata(1).mean(:,:,2);
    proj.std = std( double(greenMovie),[],3 );
    proj.max = max( greenMovie,[],3 );

    roiTemplate = struct( 'ind',[], 'pix',[], 'Npix',NaN, 'mask',[], 'edge',[], 'Nedge',NaN, 'cent',[] );
    % Start the segmentation process
    opt = {[0.06,0.03], [0.04,0.02], [0.01,0.01]};  % {[vert, horz], [bottom, top], [left, right] }
    FS = 12;
    close all; clearvars sp h;
    figure('units','normalized','WindowState','max', 'name',expt.name, 'color','k'); % 'OuterPosition',[0,0,0.48,1]
    sp(1) = subtightplot(1,3,1,opt{:});
    imshow( proj.mean, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
    title('Mean Projection', 'FontSize',FS, 'Color','w');
    sp(2) = subtightplot(1,3,2,opt{:});
    imshow( proj.max, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
    title('Max Projection', 'FontSize',FS, 'Color','w');
    sp(3) = subtightplot(1,3,3,opt{:});
    imshow( proj.std, [], 'InitialMagnification','fit' ); hold on; axis off; impixelinfo;
    title('Standard Deviation Projection', 'FontSize',FS, 'Color','w');

    % Start the segmentation process
    fprintf('\nStarting segmentation of %s', expt.name);
    % DEFINE BACKGROUND ROI
    title('Click one of the images to use for segmentation', 'FontSize',FS, 'Color','w'); waitforbuttonpress;
    title('Draw an ellipse around a background region, then press enter', 'FontSize',FS, 'Color','w');
    backEllipse = drawellipse; pause;
    backROI = roiTemplate;
    backROI.mask = backEllipse.createMask;
    backProps = regionprops( backROI.mask, 'centroid', 'pixellist', 'pixelidxlist', 'area' );
    backROI.ind = backProps.PixelIdxList;
    backROI.pix = backProps.PixelList;
    backROI.Npix = backProps.Area;
    backROI.cent = backProps.Centroid;
    backROI.edge = bwboundaries( backROI.mask, roiParam.conn, 'noholes' );
    backROI.edge = flip( backROI.edge{1} , 2);
    backROI.Nedge = numel(backROI.edge);
    plot( backROI.edge(:,1), backROI.edge(:,2), '--', 'Color', 'w' );  % outline background area

    % Threshold to separate putative signal and background
    backThresh = ceil( mean(proj.mean(backROI.ind))+roiParam.Nsigma*std( double(proj.mean(backROI.ind))) );
    tempQuant = imquantize( proj.mean, backThresh ); % imshow( tempQuant, [] );
    tempBin = logical(tempQuant - 1); % ignore the subthreshold pixels   % imshow(tempBin );
    CC = bwconncomp(tempBin,roiParam.conn);
    RP = regionprops( CC, 'Area' );
    tempThresh = CC.PixelIdxList;
    tempThresh = tempThresh( [RP.Area] > minArea ); % only keep sufficiently large ROIs
    NthreshROI = numel( tempThresh );
    allAbove = false( movieRes(1), movieRes(2) );
    threshROI = repmat(roiTemplate, 0, 1);
    for c = flip(1:NthreshROI)
        threshROI(c).ind = tempThresh{c};
        [threshROI(c).pix(:,2), threshROI(c).pix(:,1)]  = ind2sub( movieRes, tempThresh{c} );
        threshROI(c).Npix = numel(tempThresh{c});
        threshROI(c).mask = false( movieRes(1), movieRes(2) ); threshROI(c).mask( tempThresh{c} ) = true;
        threshROI(c).edge = bwboundaries( threshROI(c).mask, roiParam.conn, 'holes' );
        threshROI(c).Nedge = numel( threshROI(c).edge );
        threshROI(c).edge = cellfun( @flip, threshROI(c).edge, repmat({2},threshROI(c).Nedge,1), 'UniformOutput', false); % bwboundaries output has [y,x] columns
    end
    allAbove( vertcat( threshROI.ind ) ) = true;
    for c = 1:NthreshROI
        for n = 1:threshROI(c).Nedge
            plot( threshROI(c).edge{n}(:,1), threshROI(c).edge{n}(:,2), 'b--' );
        end
    end
    linkaxes(sp,'xy');
    % IDENTIFY CELLS
    cellROI = repmat(roiTemplate, 0, 1); somaROI = repmat(roiTemplate, 0, 1); procROI = repmat(roiTemplate, 0, 1); % initialize
    for c = 1:50 % NthreshROI % while ishandle( proj.mean ) && c <= NthreshROI %
        title('Zoom in on the cell, then unpause', 'FontSize',FS, 'Color','w'); pause;
        title(sprintf('c = %d: Press any button to continue scoring OR close the figure to stop',c), 'FontSize',FS, 'Color','w');
        try
            waitforbuttonpress;
        catch
            c = c - 1; %#ok<FXSET>
            fprintf('\n  Figure closed. Found %d good cells', c );
            break;
        end
        % Entire cell
        title(sprintf('c = %d: Draw a polygon around the entire cell territory, then unpause',c), 'FontSize',FS, 'Color','w');
        tempCell = drawpolygon; %impoly;
        pause;
        cellROI(c).mask = tempCell.createMask;%imshow(tempCell.createMask);
        cellROI(c).mask = cellROI(c).mask & allAbove; %
        cellProps = regionprops( cellROI(c).mask, 'centroid', 'pixellist', 'pixelidxlist', 'area' );
        cellROI(c).ind = cellProps.PixelIdxList;
        cellROI(c).pix = cellProps.PixelList;
        cellROI(c).Npix = cellProps.Area;
        cellROI(c).cent = cellProps.Centroid;
        cellROI(c).edge = bwboundaries( cellROI(c).mask, roiParam.conn, 'noholes' );
        cellROI(c).edge = cellfun( @flip, cellROI(c).edge , repmat({2},numel(cellROI(c).edge),1), 'UniformOutput', false);
        cellROI(c).Nedge = numel(cellROI(c).edge);
        for n = 1:cellROI(c).Nedge
            plot( cellROI(c).edge{n}(:,1), cellROI(c).edge{n}(:,2), '--', 'Color', 'c' );  % colorblind(c,:)
        end
        % Soma only
        title(sprintf('c = %d: Draw an ellipse around the soma, then unpause',c), 'FontSize',FS, 'Color','w');
        tempSoma = drawellipse; %imellipse;
        pause;
        somaROI(c).mask = tempSoma.createMask;
        somaROI(c).mask = somaROI(c).mask & allAbove; %imshow(somaMask{c});
        somaProps = regionprops( somaROI(c).mask, 'centroid', 'pixellist', 'pixelidxlist', 'area' );
        somaROI(c).ind = somaProps.PixelIdxList;
        somaROI(c).pix = somaProps.PixelList;
        somaROI(c).Npix = somaProps.Area;
        somaROI(c).cent = somaProps.Centroid;
        somaROI(c).edge = bwboundaries( somaROI(c).mask, roiParam.conn, 'noholes' );
        somaROI(c).edge = cellfun( @flip, somaROI(c).edge , repmat({2},numel(somaROI(c).edge),1), 'UniformOutput', false);
        somaROI(c).Nedge = numel(somaROI(c).edge);
        for n = 1:somaROI(c).Nedge
            plot( somaROI(c).edge{n}(:,1), somaROI(c).edge{n}(:,2), '--', 'Color', 'g' );  % colorblind(c,:)
        end
        % Processes
        procROI(c).mask = false(movieRes);
        procROI(c).mask(cellROI(c).ind) = true;
        procROI(c).mask(somaROI(c).ind) = false; %imshow(procROI(c).mask);
        procProps = regionprops( procROI(c).mask, 'centroid', 'pixellist', 'pixelidxlist', 'area' );
        if ~isempty(procProps)
            procROI(c).ind = procProps.PixelIdxList;
            procROI(c).pix = procProps.PixelList;
            procROI(c).Npix = procProps.Area;
            procROI(c).cent = procProps.Centroid;
            if numel(procROI(c).Npix) > 10
                procROI(c).edge = bwboundaries( procROI(c).mask, roiParam.conn, 'noholes' );
                procROI(c).edge  = cellfun( @flip, procROI(c).edge , repmat({2},numel(procROI(c).edge ),1), 'UniformOutput', false);
                procROI(c).Nedge = numel(procROI(c).edge);
            else
                procROI(c).edge = [NaN, NaN]; procROI(c).Nedge = NaN;
            end
        else
            procROI(c).ind = [];
            procROI(c).pix = nan(0,2);
            procROI(c).Npix = 0;
            procROI(c).cent = nan(0,2);
        end

    end
    Ncell = numel(cellROI);

    % Write cells to temporary variables for saving
    savePath = sprintf('%s%s_ROI.mat', expt.dir, expt.name); % metadata(1).root
    save(savePath, 'backROI', 'cellROI', 'somaROI', 'procROI', 'Ncell', 'proj', 'roiParam' ); fprintf('\n  Saved %s\n',savePath);
    roi2aqua( expt );
end
end