function [cellEvt, aquaOpts, aquaResults] = GetAquaResults( expt )
% Gather results of AQuA anaylysis (individual movie or full FOV) and repackage by cell and movie
IP = inputParser;
addRequired( IP, 'dataDir', @isstruct) %dataDir = 'D:\2photon\AB18\190905\008\';
parse( IP, expt ); 

% Determine whether to get results from a single movie, or all movies from an FOV
if dataDir(end) ~= filesep, dataDir = [dataDir, filesep]; end
if ~isempty(fov) 
    [~,~,tempOther] = LoadProcessed( dataDir, 'pool' ); 
    movies = tempOther(1).fovMovies;  
elseif ~isempty(movieNum)
    movies = movieNum;  
else
    error('Invalid data directory!' );
end
Nmovies = numel(movies);

cellEvt = cell(Nmovies,1); aquaOpts = cell(Nmovies,1); aquaResults = cell(Nmovies,1);
for m = flip(1:Nmovies)
    [~,movieDir] = SetDataPath( mouse, exptDate, movies(m), [] );
    aquaDir = FileFind( movieDir ); 
    if ~isempty(aquaDir)
        aquaPath = FileFind( aquaDir{1,2}, 'mat', false, @(x)(contains( x, 'AQuA' )) ); aquaPath = aquaPath{2};
        if ~isempty(aquaPath)
            [cellEvt{m}, aquaOpts{m}, aquaResults{m}] = UnpackAquaData( aquaPath );
        else
            fprintf('\nNo AQuA results file found in %s\n', aquaPath );
        end
    end
end
cellEvt = cat(1, cellEvt{:});
aquaOpts = cat(1, aquaOpts{:});
aquaResults = cat(1, aquaResults{:});

%{
figure;
for m = 1
    for c = 1
        for v = 1:cellEvt(m,c).N
            subplot(1,2,1); cla;
            imshow( cellEvt(m,c).foot(:,:,v), [] ); impixelinfo;
            title( sprintf('[c,v] = [%i,  %i]', c, v) );
            
            subplot(1,2,2);
            plot( cellEvt(m,c).dFFfull(v,:), 'r' ); hold on;
            plot( cellEvt(m,c).frames(1,v):cellEvt(m,c).frames(2,v),  cellEvt(m,c).dFF{v}, 'k' ); 
            %title( sprintf('p = %2.5f', cellEvt(m,c).dFFmaxPval(v)) );
            ylabel('dF/F'); xlabel('Frames');
            pause;
            cla;
        end
    end
end
%}       


%{
if all(cellfun(@isempty, aquaOpts)) %~exist('cellEvt','var')
    fprintf('\n%s: No AQuA results found\n', dataDir);
    cellEvt = []; %struct();
    aquaResults = []; % struct();
end
%}

end