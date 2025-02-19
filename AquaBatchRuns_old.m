function [cellEvt, aquaOpts, aquaResults] = AquaBatchRuns(expt, projParam) % , tempRes, spatRes, excBorder
%% setup
% 'projDir' is the folder containing tifs you want to deal with.
% Suggest sort the files in order, so that you can set the parameters 
% conviniently. AQuA/cfg/parameters_for_batch is the parameters excel.
% The script will read the parameters from that excel to deal with data.
% How many files you have, how many parameter settings should be in that excel.

%startup;  % initialize
%load('random_Seed.mat');
%rng(s);

% Load cell boundaries and  and landmarks/soma
[~,p_cell] = FileFinder(expt.dir, 'type','mat', 'contains','Cell'); % cell boundary path, if you have
if ~isempty(p_cell), p_cell = p_cell{1}; end
[~,p_landmark] = FileFinder(expt.dir, 'type','mat', 'contains','Soma'); % cell boundary path, if you have
if ~isempty(p_landmark), p_landmark = p_landmark{1}; end
bd = containers.Map;
bd('None') = [];
if(~strcmp(p_cell,''))
    cell_region = load(p_cell);
    bd('cell') = cell_region.bd0;
end
if(~strcmp(p_landmark,''))
    landmark = load(p_landmark);
    bd('landmk') = landmark.bd0;
end

% Get the projection parameters (projParam) and use them for Aqua options
opts = util.parseParam(1); % preset 1: original Lck. 2: Jen Lck. 5: GluSnFR
opts.frameRate = 1/projParam.rate_bin; % seconds per frame, not other way around
opts.spatialRes = projParam.umPerPixel_scale; % microns per pixel
opts.regMaskGap = 0; % how many pixels, around edge, to be excluded

% Find the projections from the green channel of individual runs
isReg = @(x)(contains(x,'_reg_'));
catProjDir = [expt.dir,'Projections']; % tif folder
[catProjNames,catProjTifs] = FileFinder(catProjDir, 'type','tif', 'contains','green', 'criteria',isReg, 'keepExt',true);
runProjDir = [expt.dir,'Projections\Runs']; % tif folder
[runProjNames,runProjTifs] = FileFinder(runProjDir, 'type','tif', 'contains','green', 'criteria',isReg, 'keepExt',true);

% Preprocess the concatenated data 
[catDatOrg,opts] = burst.prep1(catProjDir,catProjNames{1},[],opts);  % read data

% Setup the output folder
outputDir = [expt.dir,'AQuA\']; mkdir(outputDir);
cellEvt = cell(expt.Nruns,1); aquaOpts = cell(expt.Nruns,1); aquaResults = cell(expt.Nruns,1);
for runs = 1:expt.Nruns
    runName = sprintf('%s_run%i', expt.name, runs);
    fprintf('\nApplying AQuA to %s: ', runName)
    runResultPath = sprintf('%s%s_AQuA.mat', outputDir, runName);
    if ~exist(runResultPath, 'file')
        tic;
        % Load the run's projection
        %projTifName = runProjNames{contains(runProjTifs, sprintf('run%i_', runs))};
        % Note: Setting the parameters should be consistent with your target file
        %[runDatOrg,opts] = burst.prep1(runProjDir,projTifName,[],opts);  % read data
        % Use the concatenated version instead
        runDatOrg = catDatOrg(:,:,projParam.binLims(runs)+1:projParam.binLims(runs+1));
        runOpts = opts;
        runOpts.sz(3) = projParam.Nbin(runs);

        %% detection
        %sz = runOpts.sz;
        evtSpatialMask = ones(runOpts.sz(1),runOpts.sz(2));
        if bd.isKey('cell')
            bd0 = bd('cell');
            evtSpatialMask = zeros(runOpts.sz(1),runOpts.sz(2));
            for ii=1:numel(bd0)
                idx = bd0{ii}{2};
                spaMsk0 = zeros(runOpts.sz(1),runOpts.sz(2));
                spaMsk0(idx) = 1;
                evtSpatialMask(spaMsk0>0) = 1;
            end
        end
        [dat,dF,~,lmLoc,runOpts,dL] = burst.actTop(runDatOrg,runOpts,evtSpatialMask);  % foreground and seed detection
        [svLst,~,riseX] = burst.spTop(dat,dF,lmLoc,evtSpatialMask,runOpts);  % super voxel detection
        [riseLst,datR,evtLst,seLst] = burst.evtTop(dat,dF,svLst,riseX,runOpts,[],bd);  % events
        [ftsLst,dffMat] = fea.getFeatureQuick(dat,evtLst,runOpts);

        % filter by significance level
        mskx = ftsLst.curve.dffMaxZ>runOpts.zThr;
        dffMatFilterZ = dffMat(mskx,:);
        evtLstFilterZ = evtLst(mskx);
        tBeginFilterZ = ftsLst.curve.tBegin(mskx);
        riseLstFilterZ = riseLst(mskx);

        % merging (glutamate)
        if runOpts.ignoreMerge==0
            evtLstMerge = burst.mergeEvt(evtLstFilterZ,dffMatFilterZ,tBeginFilterZ,runOpts,bd);
        else
            evtLstMerge = evtLstFilterZ;
        end

        % reconstruction (glutamate)
        if runOpts.extendSV==0 || runOpts.ignoreMerge==0 || runOpts.extendEvtRe>0
            [riseLstE,datRE,evtLstE] = burst.evtTopEx(dat,dF,evtLstMerge,runOpts);
        else
            riseLstE = riseLstFilterZ; datRE = datR; evtLstE = evtLstFilterZ;
        end

        % feature extraction
        [ftsLstE,dffMatE,dMatE] = fea.getFeaturesTop(runDatOrg,evtLstE,runOpts);
        ftsLstE = fea.getFeaturesPropTop(dat,datRE,evtLstE,ftsLstE,runOpts);

        % update network features
        %sz = size(runDatOrg);
        evtx1 = evtLstE;
        ftsLstE.networkAll = [];
        ftsLstE.network = [];
        if bd.isKey('cell')
            bd0 = bd('cell');
            evtSpatialMask = zeros(runOpts.sz(1),runOpts.sz(2));
            regLst = cell(numel(bd0),1);
            for ii=1:numel(bd0)
                pix00 = bd0{ii}{2};
                regLst{ii} = pix00;
                evtSpatialMask(pix00) = 1;
            end
        else
            regLst = [];
            evtSpatialMask = ones(runOpts.sz(1),runOpts.sz(2));
        end

        if bd.isKey('landmk')
            bd1 = bd('landmk');
            lmkLst = cell(numel(bd1),1);
            for ii=1:numel(bd1)
                lmkLst{ii} = bd1{ii}{2};
            end
        else
            lmkLst = [];
        end

        try
            if ~isempty(regLst) || ~isempty(lmkLst)
                fprintf('Updating region and landmark features ...\n')
                ftsLstE.region = fea.getDistRegionBorderMIMO(evtLstE,datR,regLst,lmkLst,runOpts.spatialRes,runOpts.minShow1);
                if bd.isKey('cell')
                    bd0 = bd('cell');
                    for i = 1:numel(regLst)
                        cname{i} = bd0{i}{4};
                        if(strcmp(cname{i},'None'))
                            cname{i} = num2str(i);
                        end
                    end
                    ftsLstE.region.cell.name = cname;
                end
                if bd.isKey('landmk')
                    bd0 = bd('landmk');
                    for i = 1:numel(lmkLst)
                        lname{i} = bd0{i}{4};
                        if(strcmp(lname{i},'None'))
                            lname{i} = num2str(i);
                        end
                    end
                    ftsLstE.region.landMark.name = lname;
                end
            end
        catch
        end

        try
            ftsLstE.networkAll = fea.getEvtNetworkFeatures(evtLstE,runOpts.sz);  % all filtered events
            ftsLstE.network = fea.getEvtNetworkFeatures(evtx1,runOpts.sz);  % events inside cells only
        catch
        end


        %% export table
        fts = ftsLstE;
        tb = readtable('userFeatures.csv','Delimiter',',');
        if(isempty(ftsLstE.basic))
            nEvt = 0;
        else
            nEvt = numel(ftsLstE.basic.area);
        end
        nFt = numel(tb.Name);
        ftsTb = nan(nFt,nEvt);
        ftsName = cell(nFt,1);
        ftsCnt = 1;
        dixx = ftsLstE.notes.propDirectionOrder;
        lmkLst = [];

        for ii=1:nFt
            cmdSel0 = tb.Script{ii};
            ftsName0 = tb.Name{ii};
            % if find landmark or direction
            if ~isempty(strfind(cmdSel0,'xxLmk')) %#ok<STREMP>
                for xxLmk=1:numel(lmkLst)
                    try
                        eval([cmdSel0,';']);
                    catch
                        fprintf('Feature "%s" not used\n',ftsName0)
                        x = nan(nEvt,1);
                    end
                    ftsTb(ftsCnt,:) = reshape(x,1,[]);
                    ftsName1 = [ftsName0,' - landmark ',num2str(xxLmk)];
                    ftsName{ftsCnt} = ftsName1;
                    ftsCnt = ftsCnt + 1;
                end
            elseif ~isempty(strfind(cmdSel0,'xxDi')) %#ok<STREMP>
                for xxDi=1:4
                    try
                        eval([cmdSel0,';']);
                        ftsTb(ftsCnt,:) = reshape(x,1,[]);
                    catch
                        fprintf('Feature "%s" not used\n',ftsName0)
                        ftsTb(ftsCnt,:) = nan;
                    end
                    ftsName1 = [ftsName0,' - ',dixx{xxDi}];
                    ftsName{ftsCnt} = ftsName1;
                    ftsCnt = ftsCnt + 1;
                end
            else
                try
                    eval([cmdSel0,';']);
                    ftsTb(ftsCnt,:) = reshape(x,1,[]);
                catch
                    fprintf('Feature "%s" not used\n',ftsName0)
                    ftsTb(ftsCnt,:) = nan;
                end
                ftsName{ftsCnt} = ftsName0;
                ftsCnt = ftsCnt + 1;
            end
        end
        featureTable = table(ftsTb,'RowNames',ftsName);
        ftb = [outputDir,runName,'_FeatureTable.xlsx'];
        writetable(featureTable,ftb,'WriteVariableNames',0,'WriteRowNames',1);

        %% export movie
        datL = zeros(runOpts.sz);
        for i = 1:numel(evtLstE)
            datL(evtLstE{i}) = i;
        end
        ov1 = zeros(runOpts.sz(1),runOpts.sz(2),3,runOpts.sz(3));
        % re-scale movie
        c0 = zeros(nEvt,3);
        for nn=1:nEvt
            x = rand(1,3);
            while (x(1)>0.8 && x(2)>0.8 && x(3)>0.8) || sum(x)<1
                x = rand(1,3);
            end
            x = x/max(x);
            c0(nn,:) = x;
        end

        for tt=1:runOpts.sz(3)
            if mod(tt,100)==0
                fprintf('Frame %d\n',tt);
            end
            dat0 = runDatOrg(:,:,tt);
            if runOpts.usePG==1
                dat0 = dat0.^2;
            end
            datx = cat(3,dat0,dat0,dat0);
            datxCol = datx;
            [H,W,~] = size(datx);
            reCon = double(datRE(:,:,tt))/255;
            rPlane = zeros(H,W);
            gPlane = rPlane;
            bPlane = rPlane;
            map = datL(:,:,tt);
            rPlane(map>0) = c0(map(map>0),1);
            gPlane(map>0) = c0(map(map>0),2);
            bPlane(map>0) = c0(map(map>0),3);
            datxCol(:,:,1) = rPlane.*reCon + datxCol(:,:,1);
            datxCol(:,:,2) = gPlane.*reCon + datxCol(:,:,2);
            datxCol(:,:,3) = bPlane.*reCon + datxCol(:,:,3);
            ov1(:,:,:,tt) = datxCol;
        end
        fmov = [outputDir,runName,'_Movie.tif'];
        io.writeTiffSeq(fmov,ov1,8);


        %% export to GUI
        res = fea.gatherRes(runDatOrg,runOpts,evtLstE,ftsLstE,dffMatE,dMatE,riseLstE,datRE);
        % aqua_gui(res);
        % visualize the results in each step
        %{
        if 0
            ov1 = plt.regionMapWithData(arLst,runDatOrg,0.5); zzshow(ov1);
            ov1 = plt.regionMapWithData(svLst,runDatOrg,0.5); zzshow(ov1);
            ov1 = plt.regionMapWithData(seLst,runDatOrg,0.5,datR); zzshow(ov1);
            ov1 = plt.regionMapWithData(evtLst,runDatOrg,0.5,datR); zzshow(ov1);
            ov1 = plt.regionMapWithData(evtLstFilterZ,runDatOrg,0.5,datR); zzshow(ov1);
            ov1 = plt.regionMapWithData(evtLstMerge,runDatOrg,0.5,datR); zzshow(ov1);
            [ov1,lblMapS] = plt.regionMapWithData(evtLstE,runDatOrg,0.5,datRE); zzshow(ov1);
        end
        %}

        %% save output
        res.bd = bd;
        fprintf('Writing %s', runResultPath)
        save(runResultPath, 'res'); % path0 -> outputDir
        toc
    else
        fprintf('Already done')
    end

    % Load AQuA results
    [cellEvt{runs}, aquaOpts{runs}, aquaResults{runs}] = UnpackAquaData( runResultPath );
end

cellEvt = cat(1, cellEvt{:});
aquaOpts = cat(1, aquaOpts{:});
aquaResults = cat(1, aquaResults{:});
end