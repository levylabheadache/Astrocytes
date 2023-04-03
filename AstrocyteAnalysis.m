clear; clc; close all;
dataDir = 'D:\2photon\';  %'C:\2photon';
dataSet = 'Astrocyte'; % 'Pollen'; %   'NGC'; % 'Neutrophil'; % 
% Parse data table
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; % 'D:\MATLAB\NGCdata.xlsx'; 
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), ...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'csd',find(contains(colNames, 'CSD')), 'ref',find(contains(colNames, 'Ref')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Projection parameters
projParam.rate = 0.5; % Hz
projParam.umPerPixel = 1; % um/pix
projParam.edge = [100,100,30,30]; % crop this many pixels from each side: [left, right, top, bottom]
projParam.color = {'red','green'};
projParam.z = 1; %{3:8, 7:11};%8:10;
projParam.overwrite = true;

% Initialize variables
expt = cell(1,Nexpt); runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); 
cellEvt = cell(1,Nexpt); aquaOpts = cell(1,Nexpt); aquaResults = cell(1,Nexpt); 
Trel = cell(1,Nexpt); Fastro = cell(1,Nexpt); Fback = cell(1,Nexpt);
% Choose which subset to  process
xPresent = 55; %6; 
Npresent = numel(xPresent);
for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Parse data table, get metadata (catInfo and expt)  and Tscan
    [expt{x}, runInfo{x}] = ParseDataTable(dataTable, x, dataCol, dataDir );  
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x}); 
    [catInfo(x), expt{x}] = ConcatenateRunInfo(expt{x}, runInfo{x}, 'suffix','sbxcat'); % Get concatenated metadata
    
    % Draw the ROIs
    [backROI{x}, cellROI{x}, somaROI{x}, procROI{x}, Ncell{x}, proj{x}] = MakeAstroROI( expt{x} );

    % Run AQuA on the green projections data
    projParam = GenerateExptProjections(expt{x}, catInfo(x), Tscan{x});
    [cellEvt{x}, aquaOpts{x}, aquaResults{x}] = AquaBatchRuns(expt{x}, projParam);
    [Fastro{x}, Fback{x}, Trel{x}] = GetAstroFluor( expt{x}, runInfo{x}, projParam, backROI{x}, cellROI{x}, somaROI{x}, procROI{x}, 'overwrite',true, 'show',true );
end