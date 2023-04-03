function exptData = CombineFOVdata(Tmin, Fastro, Ncell, varargin) % velocity, 
% Combine all data from different FOV of an experiment an into a super matrix: columns [T, Fastro1, Fastro2, ... FastroN, velocity, locoVec]
IP = inputParser;
addRequired( IP, 'Tmin', @iscell )
addRequired( IP, 'Fastro', @iscell )
addRequired( IP, 'Ncell', @isnumeric )
addOptional( IP, 'velocity', {}, @iscell )
addOptional( IP, 'bout', {}, @iscell )
%addOptional( IP, 'cellEvt', {}, @iscell )
parse( IP, Tmin, Fastro, Ncell, varargin{:} ); % , velocity
velocity = IP.Results.velocity;
locoBout = IP.Results.bout;
%cellEvt = IP.Results.cellEvt;

% Determine the size of the supermatrix
NcellCum = [0, cumsum(Ncell)];
Nfov = numel( Tmin );
Nsamp = 0;
fovCol = cell(1,Nfov);
for f = 1:Nfov
    Nsamp = Nsamp + sum(cellfun( @numel, Tmin{f} ));
    fovCol{f} = NcellCum(f)+1:NcellCum(f+1);
end  
%
exptData = nan( Nsamp, sum(Ncell)+3 ); 
z = 0;
for f = 1:Nfov
    NframeMov = cellfun( @numel, Tmin{f} );
    NframeFOV = sum(NframeMov);
    exptFrames = z + (1:NframeFOV)';
    exptData(exptFrames,1) = vertcat(Tmin{f}{:}); %fovData(:,1); % Time
    exptData(exptFrames,fovCol{f}+1) = vertcat( Fastro{f}(:).cell ); % fovData(:,2:end); % fluor
    % velocity
    if ~isempty(velocity) && ~any( cellfun(@isempty, velocity{f}) ) %
        exptData(exptFrames,end-1) = vertcat(velocity{f}{:}); 
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
        exptData(exptFrames,end) = locoVec;
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
exptData = sortrows(exptData);
%plot( exptData(:,1), exptData(:,2:end-2) )
%plot( exptData(:,1), exptData(:,end-1) )

end