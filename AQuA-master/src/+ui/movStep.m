function datxCol = movStep(f,n,ovOnly,updtAll)
    % use btSt.sbs, btSt.leftView and btSt.rightView to determine what to show
    btSt = getappdata(f,'btSt');
    fh = guidata(f);
        
    if ~isfield(btSt,'GaussFilter') ||(btSt.GaussFilter==0) 
        dat = getappdata(f,'datOrg');
    else
        dat = getappdata(f,'dat');  
    end
    
    dF = getappdata(f,'dF');
    
    if isempty(dF)
        dF = dat;
    end
    
    if isempty(dat)
        dat = getappdata(f,'dat');
    end
    scl = getappdata(f,'scl');
    btSt = getappdata(f,'btSt');
    
    if ~isfield(fh,'maxPro')
        fh.maxPro = max(dat,[],3);
    end
    
    if ~exist('ovOnly','var') || isempty(ovOnly)
        ovOnly = 0;
    end
    if ~exist('n','var') || isempty(n)
        n = fh.sldMov.Value;
    end
    if ~exist('updtAll','var') || isempty(updtAll)
        updtAll = 0;
    end
    
    % re-scale movie
    dat0 = dat(:,:,n);
    dF0 = dF(:,:,n)/max(dF(:));
    if scl.map==1
        dat0 = dat0.^2;
    end
    dat0 = (dat0-scl.min)/max(scl.max-scl.min,0.01)*scl.bri;
    datx = cat(3,dat0,dat0,dat0);
    datxL = datx/scl.bri*scl.briL;
    datxR = datx/scl.bri*scl.briR;
    datxCol = ui.over.addOv(f,datx,n);
    datxColL = ui.over.addOv(f,datxL,n);
    datxColR = ui.over.addOv(f,datxR,n);
    
    if ovOnly>0
        return
    end
    
    %% overlay
    if btSt.sbs==0
        fh.ims.im1.CData = flipud(datxCol);
        ui.mov.addPatchLineText(f,fh.mov,n,updtAll);
    end
    if btSt.sbs==1
        viewName = {'leftView','rightView'};
        imName = {'im2a','im2b'};
        axLst = {fh.movL,fh.movR};
        for ii=1:2
            curType = btSt.(viewName{ii});
            axNow = axLst{ii};
            
            % clean all patches
            if updtAll>0
                types = {'quiver','line','patch','text'};
                for jj=1:numel(types)
                    h00 = findobj(axNow,'Type',types{jj});
                    if ~isempty(h00)
                        delete(h00);
                    end
                end
            end
            switch curType
                case 'Raw'
                    if ii==1
                        fh.ims.(imName{ii}).CData = flipud(datxL);
                    else
                        fh.ims.(imName{ii}).CData = flipud(datxR);
                    end
%                     ui.mov.addPatchLineText(f,axNow,n,updtAll);
                case 'dF'
                        dF0 = (dF0-scl.min)/max(scl.max-scl.min,0.01)*scl.bri;
                        dFx = cat(3,dF0,dF0,dF0);
                        dFxL = dFx/scl.bri*scl.briL;
                        dFxR = dFx/scl.bri*scl.briR;
                        if ii==1
                            fh.ims.(imName{ii}).CData = flipud(dFxL);
                        else
                            fh.ims.(imName{ii}).CData = flipud(dFxR);
                        end
                case 'Raw + overlay'
                    if ii==1
                        fh.ims.(imName{ii}).CData = flipud(datxColL);
                    else
                        fh.ims.(imName{ii}).CData = flipud(datxColR);
                    end
                    ui.mov.addPatchLineText(f,axNow,n,updtAll);
                case 'Rising map'
                    ui.mov.showRisingMap(f,imName{ii},n);
                case 'Maximum Projection'
                    if scl.map==1
                        datM = fh.maxPro.^2;
                    end
                    if ii==1
                        datML = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briL;
                        datMxL = cat(3,datML,datML,datML);
                        fh.ims.(imName{ii}).CData = flipud(datMxL);
                    else
                       datMR = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briR;
                        datMxR = cat(3,datMR,datMR,datMR);
                        fh.ims.(imName{ii}).CData = flipud(datMxR);
                    end
                    ui.mov.addPatchLineText(f,axNow,n,updtAll);
                case 'Average Projection'
                    if scl.map==1
                        datM = fh.averPro.^2;
                    end
                    if ii==1
                        datML = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briL;
                        datMxL = cat(3,datML,datML,datML);
                        fh.ims.(imName{ii}).CData = flipud(datMxL);
                    else
                       datMR = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briR;
                        datMxR = cat(3,datMR,datMR,datMR);
                        fh.ims.(imName{ii}).CData = flipud(datMxR);
                    end
                    ui.mov.addPatchLineText(f,axNow,n,updtAll);    
                
            end
        end
    end
    
    %% finish
    % adjust area to show
    if btSt.sbs==0
        fh.mov.XLim = scl.wrg;
        fh.mov.YLim = scl.hrg;
    else
        fh.movL.XLim = scl.wrg;
        fh.movL.YLim = scl.hrg;
        fh.movR.XLim = scl.wrg;
        fh.movR.YLim = scl.hrg;
    end
    
    % frame number
    ui.mov.updtMovInfo(f,n,size(dat,3));
    % f.Pointer = 'arrow';
    
end



