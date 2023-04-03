function res = saveExp(~,~,f,file0,path0,modex)
% saveExp save experiment (and export results)

fts = getappdata(f,'fts');
if ~exist('modex','var')
    modex = 0;
end
if isempty(fts)
    msgbox('Please save after event detection\n');
    return
end
if ~exist(path0,'file') && ~isempty(path0)
    mkdir(path0);    
end


%% gather results
ff = waitbar(0,'Gathering results ...');

% if do not want to detect again, do not need to save dF
vSave0 = {...  % basic variables for results analysis
    'opts','scl','btSt','ov','bd','datOrg','evt','fts','dffMat','dMat',...
    'riseLst','featureTable','userFeatures','dF',...
    };
% vSave1 = {...  % extra variables for event detection
%     'arLst','lmLoc','svLst','seLstAll','riseX','riseLstAll','evtLstAll','ftsLstAll',...
%     'dffMatAll','datRAll','evtLstFilterZ','dffMatFilterZ','tBeginFilterZ',...
%     'riseLstFilterZ','evtLstMerge','dF'...
% };
% vSave = [vSave0,vSave1];
vSave = vSave0;

res = [];
for ii=1:numel(vSave)
    v0 = vSave{ii};
    res.(v0) = getappdata(f,v0);
end

% filter features and curves
ov = getappdata(f,'ov');
ov0 = ov('Events');
xSel = ov0.sel;
btSt = getappdata(f,'btSt');
% xSel = ov0.sel;
xSelFav = false(numel(res.evt),1);
xSelFav(btSt.evtMngrMsk) = true;

res.ftsFilter = util.filterFields(fts,xSel);
res.ftsFav = util.filterFields(fts,xSelFav);
res.evtFilter = res.evt(xSel);
res.evtFav = res.evt(xSelFav);
res.dffMatFilter = res.dffMat(xSel,:,:);
res.dffMatFav = res.dffMat(xSelFav,:,:);
if ~isempty(res.dMat)
    res.dMatFilter = res.dMat(xSel,:,:);
     res.dMatFav = res.dMat(xSelFav,:,:);
end
if ~isempty(res.riseLst)  % rising map is for super events
    res.riseLstFilter = res.riseLst(xSel);
    res.riseLstFav = res.riseLst(xSelFav);
end
res.evtSelectedList = find(xSel>0);
res.evtFavList = find(xSelFav>0);

% save raw movie with 8 or 16 bits to save space
res.opts.bitNum = 16;
res.maxVal = nanmax(res.datOrg(:));  
res.datOrg = res.datOrg/res.maxVal;
dat1 = res.datOrg*(2^res.opts.bitNum-1);
res.datOrg = uint16(dat1);

res.stg.post = 1;
res.stg.detect = 1;

if modex>0
    waitbar(1,ff);
    delete(ff);
    return
end

%% export
fh = guidata(f);
opts = getappdata(f,'opts');


% btSt = getappdata(f,'btSt');
favEvtLst = btSt.evtMngrMsk;
fout = [path0,filesep,file0];
[fpath,fname,ext] = fileparts(fout);

if fh.expEvt.Value==1
    waitbar(0.25,ff,'Saving res file...');
    if isempty(ext)
        fout = [fout,'.mat'];
    end
    save(fout,'res','-v7.3');
end



% export movie
if fh.expMov.Value==1
    waitbar(0.5,ff,'Writing movie ...');
    ov1 = zeros(opts.sz(1),opts.sz(2),3,opts.sz(3));
    for tt=1:opts.sz(3)
        if mod(tt,100)==0
            fprintf('Frame %d\n',tt); 
        end
        ov1(:,:,:,tt) = ui.movStep(f,tt,1);
    end
    ui.movStep(f);
    fmov = [fpath,filesep,fname,'.tif'];
    io.writeTiffSeq(fmov,ov1,8);
end

if fh.expFt.Value==1
    waitbar(0.75,ff,'Writing feature table ...');
    % export feature table
    ftTb = getappdata(f,'featureTable');
    if isempty(ftTb)
        ui.detect.getFeatureTable(f);
        ftTb = getappdata(f,'featureTable');
    end
    cc = ftTb{:,1};

    % all selected events
    cc1 = cc(:,xSel);
    ftTb1 = table(cc1,'RowNames',ftTb.Row);
    ftb = [fpath,filesep,fname,'.csv'];
    writetable(ftTb1,ftb,'WriteVariableNames',0,'WriteRowNames',1);

    
    bd = getappdata(f,'bd');
    
    % for each region
    if ~isempty(fts.region) && isfield(fts.region.cell,'memberIdx') && ~isempty(fts.region.cell.memberIdx)
        bdcell = bd('cell');
        fpathRegion = [fpath,'\Regions'];
        if ~exist(fpathRegion,'file') && ~isempty(fpathRegion)
            mkdir(fpathRegion);    
        end

        memSel = fts.region.cell.memberIdx(xSel,:);
        for ii=1:size(memSel,2)
            mem00 = memSel(:,ii);
            Name = 'None';
            if numel(bdcell{ii})>=4
                Name = bdcell{ii}{4};
            end
            if strcmp(Name,'None')
               Name = num2str(ii); 
            end
            cc00 = cc(:,mem00>0);
            ftTb00 = table(cc00,'RowNames',ftTb.Row);
            ftb00 = [fpathRegion,filesep,fname,'_region_',Name,'.xlsx'];
            writetable(ftTb00,ftb00,'WriteVariableNames',0,'WriteRowNames',1);
        end
    end

    % for favorite events
    if ~isempty(favEvtLst)
        cc00 = cc(:,favEvtLst);
        ftTb00 = table(cc00,'RowNames',ftTb.Row);
        ftb00 = [fpath,filesep,fname,'_favorite.xlsx'];
        writetable(ftTb00,ftb00,'WriteVariableNames',0,'WriteRowNames',1);
    end

    % region and landmark map
    f00 = figure('Visible','off');
    dat = getappdata(f,'datOrg');
    dat = mean(dat,3);
    dat = dat/max(dat(:));
    Low_High = stretchlim(dat,0.001);
    dat = imadjust(dat,Low_High);
    axNow = axes(f00);
    image(axNow,'CData',flipud(dat),'CDataMapping','scaled');
    axNow.XTick = [];
    axNow.YTick = [];
    axNow.XLim = [0.5,size(dat,2)+0.5];
    axNow.YLim = [0.5,size(dat,1)+0.5];
    axNow.DataAspectRatio = [1 1 1];
    colormap gray
    ui.mov.addPatchLineText(f,axNow,0,1)
    % saveas(f00,[fpath,filesep,fname,'_landmark.fig']);
    saveas(f00,[fpath,filesep,fname,'_landmark.png'],'png');
    delete(f00);

    % rising maps
    riseLst = getappdata(f,'riseLst');
    if ~isempty(favEvtLst)
        f00 = figure('Visible','off');
        axNow = axes(f00);
        fpathRising = [fpath,filesep,'risingMaps'];
        if ~exist(fpathRising,'file')
            mkdir(fpathRising);
        end
        for ii=1:numel(favEvtLst)
            rr = riseLst{favEvtLst(ii)};
            imagesc(axNow,rr.dlyMap);
            colorbar(axNow);
            xx = axNow.XTickLabel;
            for jj=1:numel(xx)
                xx{jj} = num2str(str2double(xx{jj})+min(rr.rgw)-1);
            end
            axNow.XTickLabel = xx;
            xx = axNow.YTickLabel;
            for jj=1:numel(xx)
                xx{jj} = num2str(str2double(xx{jj})+min(rr.rgw)-1);
            end
            axNow.YTickLabel = xx;
            axNow.DataAspectRatio = [1 1 1];
            saveas(f00,[fpathRising,filesep,num2str(favEvtLst(ii)),'.png'],'png');
        end
    end
end

delete(ff);

end







