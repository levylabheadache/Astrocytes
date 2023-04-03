function roi2aqua( expt )
% Convert ROIs obtained from MakeAstroROI to AQuA format. 
% NOTE - this function assumes that 0 pixels are trimmed around the edges when data is loaded into aqua_gui.
[~,roiMat] = FileFinder( expt.dir, 'type','mat', 'contains','_ROI' );
if ~isempty( roiMat )
    load(roiMat{1}, 'cellROI', 'somaROI','Ncell' ); 
    Nz = numel(Ncell);
    % Convert from cellROI to AQuA Cell format
    bd0 = cell(1,Nz);
    for z = 1:Nz
        bd0{z} = cell(1, Ncell(z));
        for c = 1:Ncell(z)
            bd0{z}{c}{3} = 'manual';
            bd0{z}{c}{2} = cellROI{z}(c).ind; % pixels
            bd0{z}{c}{1} = {flip(cellROI{z}(c).edge{1}, 2)}; % edge
        end
        cellFileName = sprintf('%s%s_Cell.mat', expt.dir, expt.name );
        save( cellFileName, 'bd0' ); fprintf('\nSaved %s', cellFileName);
    end

    % Convert from somaROI to AQuA Landmark format
    for z = 1:Nz
        if ~isempty(somaROI{z})
            bd0{z} = cell(1, Ncell(z));
            for c = 1:Ncell(z)
                bd0{z}{c}{3} = 'manual';
                bd0{z}{c}{2} = somaROI{z}(c).ind; % pixels
                bd0{z}{c}{1} = {flip(somaROI{z}(c).edge{1}, 2)}; % edge
            end
            landFileName = sprintf('%s%s_Soma.mat', expt.dir, expt.name );
            save( landFileName, 'bd0' ); fprintf('\nSaved %s', landFileName);
        end
    end
else
    warning('%s:  No ROI file found', expt.name );
end
end