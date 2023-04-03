function flow(~,evtDat,f,op)

fh = guidata(f);
nTabTot = numel(fh.deOutTab.TabEnables);
ixTab = fh.deOutTab.Selection;
opts = getappdata(f,'opts');
skipSteps = fh.skipSteps.Value==1;
opts.skipSteps = skipSteps;
setappdata(f,'opts',opts);

% controls
if strcmp(op,'chg')
    ixTab = evtDat.NewValue;
    fh.deOutBack.Visible = 'on';
    fh.deOutNext.Visible = 'on';
    fh.deOutNext.Enable = 'on';
    fh.deOutRun.String = 'Run';
    switch ixTab
        case 1
            fh.deOutBack.Visible = 'off';
        case 4
            fh.deOutRun.String = 'Apply';
        case 5
            fh.deOutRun.String = 'Apply';
        case 7
            fh.deOutRun.String = 'Extract';
            fh.deOutNext.Visible = 'off';
    end
    if ixTab<nTabTot
        if strcmp(fh.deOutTab.TabEnables{ixTab+1},'off')
            fh.deOutNext.Enable = 'off';
        end
    end
end

% go to previous step
if strcmp(op,'back')
    if ixTab>1
        if(skipSteps && ixTab==5)
            fh.deOutTab.Selection = ixTab-4;
        else
            fh.deOutTab.Selection = ixTab-1;
        end
    end
    fh.deOutNext.Enable = 'on';
end

% run current step
if strcmp(op,'run')
    switch ixTab
        case 1
            ui.detect.actRun([],[],f);
        case 2
            ui.detect.phaseRun([],[],f);
        case 3
            ui.detect.evtRun([],[],f);
        case 4
            ui.detect.zsRun([],[],f);
        case 5
            ui.detect.mergeRun([],[],f);
        case 6
            ui.detect.evtReRun([],[],f);
        case 7
            ui.detect.feaRun([],[],f);
    end
    if ixTab<nTabTot
        fh.deOutNext.Enable = 'on';
    end
end

% go to next step
if strcmp(op,'next')
    if ixTab<nTabTot
        if(skipSteps && ixTab==1)
            fh.deOutTab.Selection = ixTab+4;
            fh.deOutTab.TabEnables{ixTab+4} = 'on';
        else
            fh.deOutTab.Selection = ixTab+1;
            fh.deOutTab.TabEnables{ixTab+1} = 'on';
        end
        
        if(fh.deOutTab.Selection==5)
            if(skipSteps)
                fh.mergeEventCorr.Visible = 'off';
                fh.mergeEventCorrText.Visible = 'off';
            else
                fh.mergeEventCorr.Visible = 'on';
                fh.mergeEventCorrText.Visible = 'on';
            end
        end
    end
end
end


