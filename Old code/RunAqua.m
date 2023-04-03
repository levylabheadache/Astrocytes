function RunAqua(dataDir, varargin)

IP = inputParser;
addRequired( IP, 'dataDir', @ischar )
addParameter( IP, 'preset', 1, @isnumeic )

parse( IP, dataDir, varargin{:} );
preset = IP.Results.preset; % 

%tifPath = 'D:\2photon\AB23\191009\011\AB23_191009_011_reg_chn2_bin1.tif';
dataDir = 'D:\2photon\AB23\191009\013\';
preset = 1;

addpath(genpath('.')); %startup;  % initialize

[ metadata, ~, ~, ~ ] = LoadProcessed( dataDir, '' );

tifName = FileFind( dataDir, 'tif', false, @(x)(contains( x, '_reg_chn2' )) );
tifName = [tifName{1},'.tif'];

opts = util.parseParam(preset,1);
opts.frameRate = 1/metadata.rate; % sec per frame?
opts.spatialRes = 1/metadata.pixPerUm;
opts.regMaskGap = 0; % how many pixels from edge of movie to crop
%{
    opts.smoXY = 1;
    opts.thrARScl = 2;
    opts.movAvgWin = 15;
    opts.minSize = 8;
    opts.regMaskGap = 0;
    opts.thrTWScl = 5;
    opts.thrExtZ = 0.5;
    opts.extendSV = 1;
    opts.cRise = 1;
    opts.cDelay = 2;
    opts.zThr = 3;
    opts.getTimeWindowExt = 10000;
    opts.seedNeib = 5;
    opts.seedRemoveNeib = 5;
    opts.thrSvSig = 1;
    opts.extendEvtRe = 0;
%}

[datOrg,opts] = burst.prep1(dataDir,tifName,[],opts);  % read data
[dat,dF,arLst,lmLoc,opts,dL] = burst.actTop(datOrg,opts);  % foreground and seed detection: ADD evtSpatialMask input from MakeAstroROI LATER
[svLst,~,riseX] = burst.spTop(dat,dF,lmLoc,[],opts);  % super voxel detection: ADD evtSpatialMask input from MakeAstroROI LATER
[riseLst,datR,evtLst,seLst] = burst.evtTop(dat,dF,svLst,riseX,opts);  % events
[ftsLst,dffMat] = fea.getFeatureQuick(datOrg,evtLst,opts);

% filter by significance level
mskx = ftsLst.curve.dffMaxZ>opts.zThr;
dffMatFilterZ = dffMat(mskx,:);
evtLstFilterZ = evtLst(mskx);
tBeginFilterZ = ftsLst.curve.tBegin(mskx);
riseLstFilterZ = riseLst(mskx);

% merging (glutamate)
evtLstMerge = burst.mergeEvt(evtLstFilterZ,dffMatFilterZ,tBeginFilterZ,opts,[]);

% reconstruction (glutamate)
if opts.extendSV==0 || opts.ignoreMerge==0 || opts.extendEvtRe>0
    [riseLstE,datRE,evtLstE] = burst.evtTopEx(dat,dF,evtLstMerge,opts);
else
    riseLstE = riseLstFilterZ; datRE = datR; evtLstE = evtLstFilterZ;
end

% feature extraction
[ftsLstE,dffMatE,dMatE] = fea.getFeaturesTop(datOrg,evtLstE,opts);
ftsLstE = fea.getFeaturesPropTop(dat,datRE,evtLstE,ftsLstE,opts);

% update network features
sz = size(datOrg);
evtx1 = evtLstE;
ftsLstE.networkAll = [];
ftsLstE.network = [];
try
    ftsLstE.networkAll = fea.getEvtNetworkFeatures(evtLstE,sz);  % all filtered events
    ftsLstE.network = fea.getEvtNetworkFeatures(evtx1,sz);  % events inside cells only
catch
end


%% export table
%fts = ftsLstE;
tb = readtable('userFeatures.csv','Delimiter',',');
nEvt = numel(ftsLstE.basic.area);
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
writetable(featureTable,ftb,'WriteVariableNames',0,'WriteRowNames',1);

end

