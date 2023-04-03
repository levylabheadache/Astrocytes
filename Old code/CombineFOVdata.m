function [Tcomb, fovComb, Fcomb, dFFcomb, Acomb, locoComb, fovCol] = CombineFOVdata(Tmin, Fastro, dFF, A, Ncell, mType, varargin) % velocity, 
% Combine all data from different FOV of an experiment an into a super matrix: columns [T, Fastro1, Fastro2, ... FastroN, velocity, locoVec]
% OUTPUT
% locoComb - combined [velocity, boutID, 
IP = inputParser;
addRequired( IP, 'Tmin', @iscell )
addRequired( IP, 'Fastro', @iscell )
addRequired( IP, 'dFF', @iscell )
addRequired( IP, 'A', @iscell )
addRequired( IP, 'Ncell', @isnumeric )
addRequired( IP, 'mType', @isstruct )
addOptional( IP, 'velocity', {}, @iscell )
addOptional( IP, 'bout', {}, @iscell )
addOptional( IP, 'periFrame', {}, @iscell )
%addOptional( IP, 'cellEvt', {}, @iscell )
addParameter( IP, 'show', false, @islogical ) % true
parse( IP, Tmin, Fastro, dFF, A, Ncell, mType, varargin{:} ); % , velocity
velocity = IP.Results.velocity;
locoBout = IP.Results.bout;
periFrame = IP.Results.periFrame;
%cellEvt = IP.Results.cellEvt;
show = IP.Results.show;

Nfov = numel( Ncell );
if show
    figure('WindowState','max');
    if isempty(velocity)
        sp(1) = subplot(2,1,1);
    else
        SP(1) = subplot(2,2,2);
        for f = 1:Nfov
            for m = 1:numel(Tmin{f})
                plot( Tmin{f}{m}, velocity{f}{m} ); hold on;
            end
        end
        sp(1) = subplot(2,2,1);
    end
    for f = 1:Nfov
        for m = 1:numel(Tmin{f})
            plot( Tmin{f}{m}, Fastro{f}(m).cell ); hold on;
        end
    end
end

% Retain only relevant movies for combination
for f = 1:Nfov
    Tmin{f} = Tmin{f}(mType(f).expt);
    Fastro{f} = Fastro{f}(mType(f).expt);
    dFF{f} = dFF{f}(mType(f).expt);
    A{f} = A{f}(mType(f).expt);
    if ~isempty(velocity), velocity{f} = velocity{f}(mType(f).expt); end
    if ~isempty(locoBout), locoBout{f} = locoBout{f}(mType(f).expt); end
end

% Determine the size of the final supermatrices
NcellCum = [0, cumsum(Ncell)];
fovCol = cell(1,Nfov);
Nsamp = 0; totCell = sum(Ncell);
for f = 1:Nfov
    Nsamp = Nsamp + sum(cellfun( @numel, Tmin{f} ));
    fovCol{f} = NcellCum(f)+1:NcellCum(f+1);
end  
%
% Combine the relevant data across FOVs
Tcomb = nan(Nsamp, 1);  fovComb = nan(Nsamp,1); 
Fcomb = struct( 'cell',nan(Nsamp, totCell), 'soma',nan(Nsamp, totCell), 'proc',nan(Nsamp, totCell) );
dFFcomb = struct( 'cell',nan(Nsamp, totCell), 'soma',nan(Nsamp, totCell), 'proc',nan(Nsamp, totCell) );
Acomb = struct( 'cell',nan(Nsamp, totCell), 'soma',nan(Nsamp, totCell), 'proc',nan(Nsamp, totCell) );
locoComb = nan(Nsamp, 3);
z = 0;
for f = 1:Nfov
    NframeMov = cellfun( @numel, Tmin{f} );
    NframeFOV = sum(NframeMov);
    exptFrames = z + (1:NframeFOV)';
    Tcomb(exptFrames,1) = vertcat(Tmin{f}{:}); 
    fovComb(exptFrames,1) = f; % vector indicating which FOV was observed at each timepoint
    Fcomb.cell(exptFrames,fovCol{f}) = vertcat( Fastro{f}(:).cell ); 
    Fcomb.soma(exptFrames,fovCol{f}) = vertcat( Fastro{f}(:).soma ); 
    Fcomb.proc(exptFrames,fovCol{f}) = vertcat( Fastro{f}(:).proc ); 
    dFFcomb.cell(exptFrames,fovCol{f}) = vertcat( dFF{f}(:).cell ); 
    dFFcomb.soma(exptFrames,fovCol{f}) = vertcat( dFF{f}(:).soma ); 
    dFFcomb.proc(exptFrames,fovCol{f}) = vertcat( dFF{f}(:).proc ); 
    Acomb.cell(exptFrames,fovCol{f}) = vertcat( A{f}(:).cell ); 
    Acomb.soma(exptFrames,fovCol{f}) = vertcat( A{f}(:).soma ); 
    Acomb.proc(exptFrames,fovCol{f}) = vertcat( A{f}(:).proc ); 
    
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
        locoVec = nan(NframeFOV, 1); runVec = nan(NframeFOV, 1); 
        for m = find( NboutMov ) %1:NmovieFOV(f)   NmovieFOV = cellfun( @numel, Tmin );
            runVec(periFrame{f}{m}.run + NframeMovCum(m)) = 1;
            for b = 1:NboutMov(m)
                boutID = str2double(sprintf('%02d', f, m, b)); % a unique identifier for each bout
                locoVec(locoBout{f}{m}(b).frames + NframeMovCum(m)) = boutID;
            end
        end
        locoComb(exptFrames,2) = locoVec;
        %locoComb(exptFrames,3) = runVec; %vertcat(periFrame{f}(:).run);
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
% Sort the supermatrices chronologically
[Tcomb, sortInd] = sortrows(Tcomb);
fovComb = fovComb(sortInd,:);
Fcomb.cell = Fcomb.cell(sortInd,:);
Fcomb.soma = Fcomb.soma(sortInd,:);
Fcomb.proc = Fcomb.proc(sortInd,:);
dFFcomb.cell = dFFcomb.cell(sortInd,:);
dFFcomb.soma = dFFcomb.soma(sortInd,:);
dFFcomb.proc = dFFcomb.proc(sortInd,:);
Acomb.cell = Acomb.cell(sortInd,:);
Acomb.soma = Acomb.soma(sortInd,:);
Acomb.proc = Acomb.proc(sortInd,:);
locoComb = locoComb(sortInd,:);

if show
    if isempty(velocity)
        sp(2) = subplot(2,1,2);
    else
        SP(2) = subplot(2,2,4);
        plot( Tcomb, locoComb(:,2) > 0 );
        ylim([-0.05,1.05]);
        linkaxes(SP,'x');
        sp(2) = subplot(2,2,3);
    end
    plot( Tcomb, Fcomb.cell );
    linkaxes(sp,'xy');
end

end