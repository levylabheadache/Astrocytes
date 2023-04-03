function [cellEvt, aquaResults] = GetAquaResults( dataDir )
% Gather results of AQuA anaylysis (individual movie or full FOV) and repackage by cell and movie
IP = inputParser;
addRequired( IP, 'dataDir', @ischar) %dataDir = 'D:\2photon\AB18\190905\008\';
parse( IP, dataDir ); 

[mouse, exptDate, movieNum, fov] = ReverseDataPath( dataDir );
if dataDir(end) ~= filesep, dataDir = [dataDir, filesep]; end
if ~isempty(fov) 
    [~,~,tempOther] = LoadProcessed( dataDir, 'pool' ); 
    movies = tempOther(1).fovMovies;
    Nmovies = numel(movies);
    for m = flip(1:Nmovies)
        [~,movieDir] = SetDataPath( mouse, exptDate, movies(m), [] );
        aquaDir = FileFind( movieDir ); 
        if ~isempty(aquaDir)
            resultPath = FileFind( aquaDir{1,2}, 'mat', false, @(x)(contains( x, 'AQuA' )) ); resultPath = resultPath{2};
            if ~isempty(resultPath)
                fprintf('Loading %s\n', resultPath );
                tempRes = load( resultPath );
                aquaResults(m) = tempRes.res;
                Ncell = size(aquaResults(m).ftsFilter.region.cell.memberIdx, 2);
                evtFrames = [aquaResults(m).ftsFilter.loc.t0; aquaResults(m).ftsFilter.loc.t1];
                for c = flip(1:Ncell)
                    cellEvt(m,c).ind = find( aquaResults(m).ftsFilter.region.cell.memberIdx(:,c) == 1 )';
                    cellEvt(m,c).N = numel( cellEvt(m,c).ind );
                    cellEvt(m,c).area = aquaResults(m).ftsFilter.basic.area( cellEvt(m,c).ind );
                    cellEvt(m,c).dur = aquaResults(m).ftsFilter.curve.width11( cellEvt(m,c).ind );
                    cellEvt(m,c).dFF = aquaResults(m).ftsFilter.curve.dffMax2( cellEvt(m,c).ind );%
                    cellEvt(m,c).frames = evtFrames(:,cellEvt(m,c).ind);
                end
            else
                fprintf('\nNo AQuA results file found in %s\n', resultPath );
            end
        end
    end
elseif ~isempty(movieNum)
    aquaDir = FileFind( dataDir );
    if ~isempty(aquaDir)
        resultPath = FileFind( aquaDir{1,2}, 'mat', false, @(x)(contains( x, 'AQuA' )) ); resultPath = resultPath{2};
        if ~isempty(resultPath)
            fprintf('\nLoading %s', resultPath );
            tempRes = load( resultPath );
            aquaResults = tempRes.res;
            Ncell = size(aquaResults.ftsFilter.region.cell.memberIdx, 2);
            evtFrames = [aquaResults.ftsFilter.loc.t0; aquaResults.ftsFilter.loc.t1];
            for c = flip(1:Ncell)
                cellEvt(1,c).ind = find( aquaResults.ftsFilter.region.cell.memberIdx(:,c) == 1 )'; % which events are associated with this cell?
                cellEvt(1,c).N = numel( cellEvt(1,c).ind );
                cellEvt(1,c).area = aquaResults.ftsFilter.basic.area( cellEvt(1,c).ind );
                cellEvt(1,c).dur = aquaResults.ftsFilter.curve.width11( cellEvt(1,c).ind );
                cellEvt(1,c).dFF = aquaResults.ftsFilter.curve.dffMax2( cellEvt(1,c).ind );%
                cellEvt(1,c).frames = evtFrames(:,cellEvt(1,c).ind);
                
                aquaResults.evt{21} % cellEvt(1,c).ind
                
            end
        else
            fprintf('\nNo AQuA results file found in %s\n', resultPath );
        end
    end
else
    error('Invalid data directory!' );
end
if ~exist('cellEvt','var')
    fprintf('\n%s: No AQuA results found\n', dataDir);
    cellEvt = []; %struct();
    aquaResults = []; % struct();
end

end