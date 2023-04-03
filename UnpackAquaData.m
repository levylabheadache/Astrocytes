function [cellEvent, aquaOpts, aquaResults, Ncell] = UnpackAquaData(aquaPath)
% Unpacks the results of interest from the complicated AQuA results structure
% INPUTS
% aquaDir = file path to aqua results .mat file
% OUTPUTS
% cellEvt = structure array summarizing event features for each cell
% aquaOpts = AQuA parameters used to analyze each movie
% aquaResults = original ("packed") AQuA results structure
%
fprintf('\nLoading %s\n', aquaPath );
tempRes = load( aquaPath );
aquaResults = tempRes.res;
aquaOpts = aquaResults.opts;
movieDim = aquaOpts.sz;
falseMat = false(movieDim);

try
    Ncell = size(aquaResults.fts.region.cell.memberIdx, 2);
catch
    %fprintf('No cells drawn');
    Ncell = 0;
end

evtFrames = [aquaResults.fts.loc.t0; aquaResults.fts.loc.t1];
if Ncell > 0
    for c = flip(1:Ncell)
        % Basic event properties
        cellEvent(1,c).ind = find( aquaResults.fts.region.cell.memberIdx(:,c) == 1 )'; % which events are associated with this cell?
        cellEvent(1,c).N = numel( cellEvent(1,c).ind );
        cellEvent(1,c).area = aquaResults.fts.basic.area( cellEvent(1,c).ind );
        cellEvent(1,c).dur = aquaResults.fts.curve.width11( cellEvent(1,c).ind ); % time, in seconds, from crossing above 10% of peak and back down
        cellEvent(1,c).dFFmax = aquaResults.fts.curve.dffMax2( cellEvent(1,c).ind ); % maximum dF/F, excluding effects of other events (see manual)
        cellEvent(1,c).dFFmaxPval = aquaResults.fts.curve.dffMaxPval( cellEvent(1,c).ind );
        cellEvent(1,c).frames = evtFrames(:,cellEvent(1,c).ind);
        cellEvent(1,c).circ = aquaResults.fts.basic.circMetric( cellEvent(1,c).ind ); % measure of the circularity of the event's footprint
        if ~isempty(aquaResults.fts.region.landmarkDist)
            cellEvent(1,c).dist2soma = aquaResults.fts.region.landmarkDist.distMin(cellEvent(1,c).ind,c); % miminum distance to soma
        else
            cellEvent(1,c).dist2soma = NaN;
        end
        % Get the dF/F traces, excluding effects of other events (see manual)
        cellEvent(1,c).dFFfull = double(aquaResults.dffMat(cellEvent(1,c).ind,:,2)); % trace for all frames from the event's pixels
        for v = 1:cellEvent(1,c).N
            cellEvent(1,c).dFF{v} = cellEvent(1,c).dFFfull(v,cellEvent(1,c).frames(1,v):cellEvent(1,c).frames(2,v)); % event-specific dF/F values
        end
        % Get the event's footprint
        cellEvent(1,c).foot = zeros( movieDim(1), movieDim(2), cellEvent(1,c).N );
        for v = 1:cellEvent(1,c).N
            %cellEvent(1,c).foot = aquaResults.fts.basic.map{cellEvent(1,c).ind(v)};
            tempMat = falseMat;
            tempMat( aquaResults.evt{cellEvent(1,c).ind(v)} ) = true;
            cellEvent(1,c).foot(:,:,v) = sum( tempMat(:,:,cellEvent(1,c).frames(1,v):cellEvent(1,c).frames(2,v)), 3);  %or(
        end
        % Propagation
    end
else
    % Basic event properties
    cellEvent.ind = 1:numel(aquaResults.evt); %find( aquaResults.fts.region.cell.memberIdx(:,c) == 1 )'; % which events are associated with this cell?
    cellEvent.N = numel(aquaResults.evt);
    cellEvent.area = aquaResults.fts.basic.area;
    cellEvent.dur = aquaResults.fts.curve.width11; % time, in seconds, from crossing above 10% of peak and back down
    cellEvent.dFFmax = aquaResults.fts.curve.dffMax2; % maximum dF/F, excluding effects of other events (see manual)
    cellEvent.dFFmaxPval = aquaResults.fts.curve.dffMaxPval;
    cellEvent.frames = evtFrames;
    cellEvent.circ = aquaResults.fts.basic.circMetric; % measure of the circularity of the event's footprint
    if isfield(aquaResults.fts, 'region') && ~isempty(aquaResults.fts.region.landmarkDist)
        cellEvent.dist2soma = aquaResults.fts.region.landmarkDist.distMin; % miminum distance to soma
    else
        cellEvent.dist2soma = NaN;
    end
    %cellEvent.dist2soma = NaN; %aquaResults.fts.region.landmarkDist.distMin(cellEvent.ind,c); % miminum distance to soma
    % Get the dF/F traces, excluding effects of other events (see manual)
    cellEvent.dFFfull = double(aquaResults.dffMat(:,:,2)); % trace for all frames from the event's pixels
    for v = 1:cellEvent.N
        cellEvent.dFF{v} = cellEvent.dFFfull(v,cellEvent.frames(1,v):cellEvent.frames(2,v)); % event-specific dF/F values
    end
    % Get the event's footprint
    cellEvent.foot = zeros( movieDim(1), movieDim(2), cellEvent.N );
    for v = 1:cellEvent.N
        %cellEvent.foot = aquaResults.fts.basic.map{cellEvent.ind(v)};
        tempMat = falseMat;
        tempMat( aquaResults.evt{cellEvent.ind(v)} ) = true;
        cellEvent.foot(:,:,v) = sum( tempMat(:,:,cellEvent.frames(1,v):cellEvent.frames(2,v)), 3);  %or(
    end
end
end