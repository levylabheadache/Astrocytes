function [Tcomb, FcellComb, FsomaComb, FprocComb, locoComb] = CombineFOVdata(Tmin, Fastro, Ncell, mExpt, varargin) % velocity, 
% Combine all data from different FOV of an experiment an into a super matrix: columns [T, Fastro1, Fastro2, ... FastroN, velocity, locoVec]
IP = inputParser;
addRequired( IP, 'Tmin', @iscell )
addRequired( IP, 'Fastro', @iscell )
addRequired( IP, 'Ncell', @isnumeric )
addRequired( IP, 'mExpt', @iscell )
addOptional( IP, 'velocity', {}, @iscell )
addOptional( IP, 'bout', {}, @iscell )
%addOptional( IP, 'cellEvt', {}, @iscell )
addParameter( IP, 'show', false, @islogical )
parse( IP, Tmin, Fastro, Ncell, mExpt, varargin{:} ); % , velocity
velocity = IP.Results.velocity;
locoBout = IP.Results.bout;
%cellEvt = IP.Results.cellEvt;
show = IP.Results.show;

% Restrict combination to relevant movies only
Nfov = numel( Tmin );
for f = 1:Nfov


% Determine the size of the supermatrix
NcellCum = [0, cumsum(Ncell)];

fovCol = cell(1,Nfov);
Nsamp = 0;
for f = 1:Nfov
    Nsamp = Nsamp + sum(cellfun( @numel, Tmin{f} ));
    fovCol{f} = NcellCum(f)+1:NcellCum(f+1);
end  
%
% Combine the relevant data across FOVs
Tcomb = nan(Nsamp, 1); 
FcellComb = nan(Nsamp, sum(Ncell)); FsomaComb = nan(Nsamp, sum(Ncell)); FprocComb = nan(Nsamp, sum(Ncell));
locoComb = nan(Nsamp, 2);
z = 0;
for f = 1:Nfov
    NframeMov = cellfun( @numel, Tmin{f} );
    NframeFOV = sum(NframeMov);
    exptFrames = z + (1:NframeFOV)';
    Tcomb(exptFrames,1) = vertcat(Tmin{f}{:}); %fovData(:,1); % Time
    FcellComb(exptFrames,fovCol{f}) = vertcat( Fastro{f}(:).cell ); 
    FsomaComb(exptFrames,fovCol{f}) = vertcat( Fastro{f}(:).soma ); 
    FprocComb(exptFrames,fovCol{f}) = vertcat( Fastro{f}(:).proc ); 
    
    % velocity
    if ~isempty(velocity) && ~any( cellfun(@isempty, velocity{f}) ) %
        locoComb(exptFrames,1) = vertcat(velocity{f}{:}); 
    else
        fprintf('\n  No velocity data')
    end
    
    % locomotive bouts
    if ~isempty(locoBout) 
        NboutMov = cellfun( @numel, locoBout{f} );
        NframeMovCum = [0,cumsum(NframeMov)];
        locoVec = zeros(NframeFOV, 1);
        for m = find( NboutMov ) %1:NmovieFOV(f)   NmovieFOV = cellfun( @numel, Tmin );
            for b = 1:NboutMov(m)
                boutID = str2double(sprintf('%02d', f, m, b)); % a unique identifier for each bout
                locoVec(locoBout{f}{m}(b).frames + NframeMovCum(m)) = boutID;
            end
        end
        locoComb(exptFrames,2) = locoVec;
    else
        fprintf('\n  No bouts')
    end
    % calcium events (AQuA)
    %{
    if ~isempty(cellEvt) 
        evtMat = 
        for c = 1:Ncell(f)
            
        end
    else
        fprintf('\n  No AQuA events')
    end
    %}
    %
    z = z + NframeFOV;
end
% Sort the matrices by time
[Tcomb, sortInd] = sortrows(Tcomb);
FcellComb = FcellComb(sortInd,:);
FsomaComb = FsomaComb(sortInd,:);
FprocComb = FprocComb(sortInd,:);
locoComb = locoComb(sortInd,:);

%plot( exptData(:,1), exptData(:,2:end-2) )
%plot( exptData(:,1), exptData(:,end-1) )

end