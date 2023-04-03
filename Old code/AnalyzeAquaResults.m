% LOAD AQUA RESULTS FOR EACH MOVIE, ORGANIZED BY FOV
clear; clc; close all;
mouse = 'AB23'; %'AB13'; %'AB17'; % 'DL159';
date = '191009';
mainDir = 'D:\2photon\'; %'D:\Andy\AstroCa\';
xlsPath = sprintf('%s%s\\%s.xlsx', mainDir, mouse, mouse );
dateDir = sprintf('%s%s\\%s\\%s\\', mainDir, mouse, date );  % [mouseFolders{dateFoldInd,2}, movieName(end-2:end),'\'];

movieTable = readtable( xlsPath, 'sheet', date); % read the mouse's Excel file
prdCol = find( strcmpi( movieTable{1,:}, 'Prd') );  %'Period'
PlaneRows = find( cellfun( @str2double,movieTable{:,prdCol} ) == 0 ); % find all 2D movies
PlaneMovies = cellfun(@str2double, movieTable{PlaneRows, 1} );
VolRows = find( cellfun( @str2double,movieTable{:,prdCol} ) > 0 ); % find all volume scans
VolMovies = cellfun(@str2double, movieTable{VolRows, 1} );

fovCol = find( strcmpi( movieTable{1,:}, 'FoV') );
fovVec = cellfun( @str2double, movieTable{:,fovCol} );
fovList = unique( fovVec )'; fovList( isnan(fovList ) ) = [];
Nfov = numel(fovList);

movieDir = cell(1,Nfov);
fovName = cell(1,Nfov); fovDir = cell(1,Nfov);
fovPlaneMovies = cell(1,Nfov); NfovPlaneMov = zeros(1,Nfov);
metadata = cell(1,Nfov); Tadj = cell(1,Nfov);
for f = 1:Nfov
    [fovName{f}, fovDir{f}] = SetDataPath( mouse, date, [], f ); %
    fovPlaneMovies{f} = intersect( str2double( movieTable{ fovVec == f , 1  } ), PlaneMovies )'; % str2double( movieTable{ fovVec == f , 1  } )';
    NfovPlaneMov(f) = numel(fovPlaneMovies{f});
    movieDir{f} = cell(NfovPlaneMov(f), 1);
    
    %[ metadata(f, Tout, otherData, loadMovie, varargout ] = LoadProcessed( fovDir{f}, 'pool' );
    [Fastro{f}, Fback{f}, Tadj{f}] = GetAstroFluor( fovDir{f}, 'show',true, 'refTime', datetime(2019,11,18,12,56,00) ); 
    for m = flip( 1:NfovPlaneMov(f) ) %2: %fovMovies
        movieDir{f}{m} = sprintf('%s%03d\\', dateDir, fovPlaneMovies{f}(m) );
        [~, T{f,m}, ~, ~] = LoadProcessed( movieDir{f}{m}, 'reg' ); % metadata(f,m)
        aquaDir = FileFind( movieDir{f}{m} ); aquaDir = aquaDir{1,2};
        resultPath = FileFind( aquaDir, 'mat', false, @(x)(contains( x, 'AQuA' )) ); resultPath = resultPath{2};
        tempRes = load( resultPath );
        aquaResults(f,m) = tempRes.res;
    end
    
end

%% CALCULATE/PLOT BASIC RESULTS
movieDur = cell(1,Nfov); % nan(Nfov, max(NfovPlaneMov));
for f = 2
    Ncell(f) = numel( aquaResults(f,1).ftsFilter.region.cell.mask );
    cellEvtInd{f} = cell(NfovPlaneMov(f), Ncell(f) ); cellEvtArea{f} = cell(NfovPlaneMov(f), Ncell(f) ); cellEvtDur{f} = cell(NfovPlaneMov(f), Ncell(f) ); 
    cellEvtMaxdFF{f} = cell(NfovPlaneMov(f), Ncell(f) );
    for m = 1:NfovPlaneMov(f)
        movieDur{f}(m) = T{f,m}(end);
        for c = 1:Ncell(f)
            cellEvtInd{f}{m,c} = find( aquaResults(f,m).ftsFilter.region.cell.memberIdx(:,c) == 1 )';
            cellEvtArea{f}{m,c} = aquaResults(f,m).ftsFilter.basic.area( cellEvtInd{f}{m,c} );
            cellEvtDur{f}{m,c} = aquaResults(f,m).ftsFilter.curve.width11( cellEvtInd{f}{m,c} );
            cellEvtMaxdFF{f}{m,c} = aquaResults(f,m).ftsFilter.curve.dffMax2( cellEvtInd{f}{m,c} );
        end
    end
    NcellEvt{f} = cellfun( @numel, cellEvtInd{f} );
    cellEvtRate{f} = 60*NcellEvt{f}./repmat( movieDur{f}', 1, Ncell(f) ); % events/min for each cell/movie
    cellEvtMeanArea{f} = cellfun( @mean, cellEvtArea{f} ); % mean event size (square um) for each cell/movie
    cellEvtMeanDur{f} = cellfun( @mean, cellEvtDur{f} ); % mean event size (square um) for each cell/movie
    cellEvtMaxdFF{f} = cellfun( @mean, cellEvtMaxdFF{f} ); % mean event peak dF/Fo for each cell/movie
end
%%
injTimestamp = datetime('09-Oct-2019 13:03:00');
close all; clearvars sp h;
figure('units','normalized','WindowState','max', 'color','w');

fovColor = distinguishable_colors(Nfov);

Tcno = minutes([metadata(f,:).timestamp] - injTimestamp);
for f = 2 %1:2
    % Event Rate
    subplot(4,1,1);
    errorbar( Tcno, nanmean(cellEvtRate{f},2), nanstd(cellEvtRate{f},[],2)/sqrt(Ncell(f)), 'ko' ); hold on;
    plot( Tcno, cellEvtRate{f} ); hold on;
    ylabel('Event Rate (events/minute)'); 
    set(gca,'XtickLabel',[]);
    % Event Size
    subplot(4,1,2);
    errorbar( Tcno, nanmean(cellEvtMeanArea{f},2), nanstd(cellEvtMeanArea{f},[],2)/sqrt(Ncell(f)), 'ko' ); hold on;
    plot( Tcno, cellEvtMeanArea{f} ); hold on;
    ylabel('Event Area (um^2)'); 
    set(gca,'Yscale','log');
    
    % Event Duration
    subplot(4,1,3);
    errorbar( Tcno, nanmean(cellEvtMeanDur{f},2), nanstd(cellEvtMeanDur{f},[],2)/sqrt(Ncell(f)), 'ko' ); hold on;
    plot( Tcno, cellEvtMeanDur{f} ); hold on;
    ylabel('Event Duration (s)'); 
    set(gca,'Yscale','log');
    
    % Event dFFmax
    subplot(4,1,4);
    errorbar( Tcno, nanmean(cellEvtMaxdFF{f},2), nanstd(cellEvtMaxdFF{f},[],2)/sqrt(Ncell(f)), 'ko' ); hold on;
    plot( Tcno, cellEvtMaxdFF{f} ); hold on;
    ylabel('dF/F'); 
    %set(gca,'Yscale','log');
end
xlabel('Time post-CNO (min)');
%xlim([-Inf,Inf]); 


%% Raster plot for each FOV/movie/cell
movieRaster = cell(1,Nexpt);  fovRaster = cell(1,Nexpt);
for x = 1
    for f = flip(1:Nfov(x))
        for m = mShow{x}{f}
            movieRaster{x}{f}{m} = zeros(size( Fastro{x}{f}(m).cell ));
            for c = 1:Ncell{x}(f)
                for v = 1:cellEvt{x}{f}(m,c).N
                    tempFrames = cellEvt{x}{f}(m,c).frames(1,v):cellEvt{x}{f}(m,c).frames(2,v);
                    movieRaster{x}{f}{m}(tempFrames, c) = movieRaster{x}{f}{m}(tempFrames, c) + 1;
                end
            end
            %imagesc(cellRaster')
        end
        fovRaster{x}{f} = vertcat( movieRaster{x}{f}{:} );
        imagesc(fovRaster{x}{f}')
    end
end


%%

cellEvt{x}{f} = GetAquaResults( dataDir{x}{f} );


