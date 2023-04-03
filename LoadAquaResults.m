function [cellEvent, aquaOpts, aquaResults, Ncell] = LoadAquaResults(expt)

aquaDir = [expt.dir, 'AQuA\'];

cellEvent = cell(1,expt.Nruns);  aquaOpts = cell(1,expt.Nruns);  aquaResults = cell(1,expt.Nruns);  Ncell = cell(1,expt.Nruns);
for runs = 1:expt.Nruns
    [~,aquaRunPath] = FileFinder(aquaDir, 'contains',sprintf('run%i',runs), 'type','mat');
    if ~isempty(aquaRunPath)
        fprintf('\nUnpacking %s', aquaRunPath{1});
        [cellEvent{runs}, aquaOpts{runs}, aquaResults{runs}, Ncell{runs}] =  UnpackAquaData(aquaRunPath{1});
    end
end
end