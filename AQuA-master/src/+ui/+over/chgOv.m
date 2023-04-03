function chgOv(~,~,f,op)
% overlay related functions

fh = guidata(f);

% read overlay features
if op==0
    tb = readtable('userFeatures.csv','Delimiter',',');
    setappdata(f,'userFeatures',tb);
    fh.overlayFeature.String = tb.Name;
    fprintf('Reading done.\n');
    return
end

% enable or disable feature overlay
if op==1
    ovName = fh.overlayDat.String{fh.overlayDat.Value};
    if strcmp(ovName,'Events')
        xxx = 'on';
    else        
        xxx = 'off';
    end
    fh.overlayFeature.Enable = xxx;
    fh.overlayColor.Enable = xxx;
    fh.overlayTrans.Enable = xxx;
    fh.overlayScale.Enable = xxx;
    fh.overlayPropDi.Enable = xxx;
    fh.overlayLmk.Enable = xxx;
    fh.sldMinOv.Enable = xxx;
    fh.sldMaxOv.Enable = xxx;
%     fh.sldBriOv.Enable = xxx;
    fh.updtFeature1.Enable = xxx;
    return
end

% calcuate overlay
idx = fh.overlayDat.Value;
ovSel = fh.overlayDat.String{idx};

btSt = getappdata(f,'btSt');
btSt.overlayDatSel = ovSel;
btSt.overlayColorSel = 'Random';

% update color code for events
if strcmp(ovSel,'Events')    
    ovFea = fh.overlayFeature.String{fh.overlayFeature.Value};
    ovCol = fh.overlayColor.String{fh.overlayColor.Value};    
    if strcmp(ovFea,'Index')
        ovCol = 'Random';
        fh.overlayColor.Value = 1;
    else
        if strcmp(ovCol,'Random')
            ovCol = 'GreenRed';
            fh.overlayColor.Value = 2;
        end
    end    
    
    btSt.overlayFeatureSel = ovFea;
    btSt.overlayColorSel = ovCol;
        
    xxTrans = fh.overlayTrans.String{fh.overlayTrans.Value};  % transform
    xxScale = fh.overlayScale.String{fh.overlayScale.Value};  % scale
    xxDi = fh.overlayPropDi.Value;  % direction
    xxLmk = str2double(fh.overlayLmk.String);  % landmark
        
    fts = getappdata(f,'fts');
    tb = getappdata(f,'userFeatures');
    xSel = cellfun(@(x) strcmp(x,ovFea), tb.Name);
    cmdSel = tb.Script{xSel};
    if isfield(fts,'locAbs')
        nEvt = numel(fts.locAbs);
    else
        nEvt = numel(fts.basic.area);
    end
    
    % change overlay value according to user input
    try
        cVal = ui.over.getVal(fts,cmdSel,xxTrans,xxScale,xxDi,xxLmk);
    catch
        msgbox('Invalid script');
        return
    end
    
    if sum(~isnan(cVal))==0
        msgbox('This feature is not used');
        return
    end
    
    if sum(~isinf(cVal))==0
        cVal = zeros(numel(cVal),1);
    end
    
    cVal(isinf(cVal) & cVal>0) = max(cVal(~isinf(cVal)));
    cVal(isinf(cVal) & cVal<0) = min(cVal(~isinf(cVal)));
    
    % update overlay color
    [col0,cMap0] = ui.over.getColorCode(f,nEvt,ovCol,cVal);
    btSt.mapNow = cMap0;
    ov = getappdata(f,'ov');
    ov0 = ov('Events');
    ov0.col = col0;
    ov0.colVal = cVal;
    ov('Events') = ov0;
    setappdata(f,'ov',ov);
    
    % update min, max and brightness slider
    scl = getappdata(f,'scl');
    scl.minOv = min(cVal);
    scl.maxOv = max(cVal);
    setappdata(f,'scl',scl);
    
    fh.sldMinOv.Min = nanmin(cVal);
    fh.sldMinOv.Max = nanmax(cVal);
    fh.sldMinOv.Value = nanmin(cVal);
    fh.sldMaxOv.Min = nanmin(cVal);
    fh.sldMaxOv.Max = nanmax(cVal);
    fh.sldMaxOv.Value = nanmax(cVal);
    fh.sldMinOv.Enable = 'on';
    fh.sldMaxOv.Enable = 'on';
    fh.txtMinOv.String = ['Min:',num2str(scl.minOv)];
    fh.txtMaxOv.String = ['Max:',num2str(scl.maxOv)];
else
    % re-shuffle color code in other cases
    if ~strcmp(ovSel,'None')
        ov = getappdata(f,'ov');
        ov0 = ov(ovSel);
        nEvt = numel(ov0.colVal);
        col0 = ui.over.getColorCode(f,nEvt,'Random');
        ov0.col = col0;
        ov(ovSel) = ov0;
        setappdata(f,'ov',ov);
    end
end

setappdata(f,'btSt',btSt);

% update color map
ui.over.adjMov([],[],f,1);

% show movie with overlay
ui.movStep(f);
end







