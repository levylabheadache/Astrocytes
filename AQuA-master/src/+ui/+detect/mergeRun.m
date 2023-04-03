function mergeRun(~,~,f)
% z scores and filtering

fprintf('Merging ...\n')

fh = guidata(f);
ff = waitbar(0,'Merging ...');
opts = getappdata(f,'opts');
bd = getappdata(f,'bd');

try
    opts.ignoreMerge = fh.ignoreMerge.Value==1;
    opts.mergeEventDiscon = str2double(fh.mergeEventDiscon.String);
    opts.mergeEventCorr = str2double(fh.mergeEventCorr.String);
    opts.mergeEventMaxTimeDif = str2double(fh.mergeEventMaxTimeDif.String);
    setappdata(f,'opts',opts);
catch
    msgbox('Error setting parameters')
end
    
if(~opts.skipSteps)
    evtLstFilterZ = getappdata(f,'evtLstFilterZ');
    dffMatFilterZ = getappdata(f,'dffMatFilterZ');
    tBeginFilterZ = getappdata(f,'tBeginFilterZ');
    
    if opts.ignoreMerge==0
        evtLstMerge = burst.mergeEvt(evtLstFilterZ,dffMatFilterZ,tBeginFilterZ,opts,bd);
    else
        evtLstMerge = evtLstFilterZ;
    end
else
    evtLstFilterZ = getappdata(f,'arLst');
    setappdata(f,'evtLstFilterZ',evtLstFilterZ);
    setappdata(f,'datRAll',uint8(ones(opts.sz)*255));
    if opts.ignoreMerge==0
        evtLstMerge = burst.mergeEvt_SkipSteps(evtLstFilterZ,opts,bd);
    else
        evtLstMerge = evtLstFilterZ;
    end
end
setappdata(f,'evtLstMerge',evtLstMerge);

waitbar(1,ff);

ui.detect.postRun([],[],f,evtLstMerge,[],'Step 5: events merged');

fh.nEvtName.String = 'nEvt';
fh.nEvt.String = num2str(numel(evtLstMerge));
delete(ff);
fprintf('Done\n')

end





