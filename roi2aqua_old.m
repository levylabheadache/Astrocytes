function roi2aqua( expt )
% Convert ROIs obtained from MakeAstroROI to AQuA format. 
% NOTE - this function assumes that 0 pixels are trimmed around the edges when data is loaded into aqua_gui.
[~,roiMat] = FileFinder( expt.dir, 'type','mat', 'contains','ROI' );
if ~isempty( roiMat )
    roiStruct = load(roiMat{1} );
    % Convert from cellROI to AQuA Cell format
    bd0 = cell(1, roiStruct.Ncell);
    for c = 1:roiStruct.Ncell
        bd0{c}{3} = 'manual';
        bd0{c}{2} = roiStruct.cellROI(c).ind; % pixels
        bd0{c}{1} = {flip(roiStruct.cellROI(c).edge{1}, 2)}; % edge
    end
    cellFileName = sprintf('%s%s_Cell.mat', expt.dir, expt.name );
    save( cellFileName, 'bd0' ); fprintf('\nSaved %s', cellFileName);
    
    % Convert from cellROI to AQuA Landmark format
    bd0 = cell(1, roiStruct.Ncell);
    for c = 1:roiStruct.Ncell
        bd0{c}{3} = 'manual';
        bd0{c}{2} = roiStruct.somaROI(c).ind; % pixels
        bd0{c}{1} = {flip(roiStruct.somaROI(c).edge{1}, 2)}; % edge
    end
    landFileName = sprintf('%s%s_Soma.mat', expt.dir, expt.name );
    save( landFileName, 'bd0' ); fprintf('\nSaved %s', landFileName);
else
    warning('%s:  No ROI file found', dataName );
end
end