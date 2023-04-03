function [cellROI, somaROI, procROI, Ncell] = DrawAstroROI(varargin)
    
IP = inputParser;

addParameter( IP, 'conn', 8, @isnumeric ) %conn = 8;
addParameter( IP, 'radius', 5, @isnumeric ) %minRad = 5; % um
addParameter( IP, 'Nsigma', 3, @isnumeric ) % Nsigma = 2.5;
parse( IP, varargin{:} );  
roiParam.conn = IP.Results.conn;
roiParam.radius = IP.Results.radius;
roiParam.Nsigma = IP.Results.Nsigma;
%setSave = IP.Results.save; 
%setLoad = IP.Results.load; 

%opt = {[0.06,0.03], [0.04,0.02], [0.01,0.01]};  % {[vert, horz], [bottom, top], [left, right] } 
FS = 12;
roiTemplate = struct( 'ind',[], 'pix',[], 'Npix',NaN, 'mask',[], 'edge',[], 'Nedge',NaN, 'cent',[] );
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
    cellROI(c).mask = cellROI(c).mask; %& allAbove; %
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
    somaROI(c).mask = somaROI(c).mask; % & allAbove; %imshow(somaMask{c}); 
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
    procROI(c).mask = false( size(cellROI(1).mask) ); 
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

end

