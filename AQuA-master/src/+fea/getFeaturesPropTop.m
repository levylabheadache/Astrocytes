function ftsLst = getFeaturesPropTop(~,evtRec,evtLst,ftsLst,opts)
% getFeaturesPropTop extract propagation related features
% dat: single (0 to 1)
% evtMap: single ( integer)
% evtRec: uint8 ( integer)

[H,W,T] = size(evtRec);

if ~isfield(opts,'northx')
    northDi = [0,1];
else
    northDi = [opts.northx,opts.northy];
end

% if opts.usePG
%     dat = dat.^2;
% end
muPix = opts.spatialRes;
ftsLst.propagation = [];
for ii=1:numel(evtLst)
    if mod(ii,100)==0
        fprintf('%d/%d\n',ii,numel(evtLst))
    end
    pix0 = evtLst{ii};
    if isempty(pix0)
        voxi = ones(4);  % dummy events
        voxr = ones(4)*0.5;
        %continue
    else
        [ih,iw,it] = ind2sub([H,W,T],pix0);
        rgH = max(min(ih)-1,1):min(max(ih)+1,H);
        rgW = max(min(iw)-1,1):min(max(iw)+1,W);
        rgT = ftsLst.curve.tBegin(ii):ftsLst.curve.tEnd(ii);
        
        % basic and propagation features
        ih1 = ih-min(rgH)+1;
        iw1 = iw-min(rgW)+1;
        it1 = it-min(rgT)+1;
%         voxd = dat(rgH,rgW,rgT);
        voxi = zeros(length(rgH),length(rgW),length(rgT));
        pix1 = sub2ind(size(voxi),ih1,iw1,it1);
        voxi(pix1) = 1;
        voxr = evtRec(rgH,rgW,rgT);
        voxr = double(voxr)/255;
    end
    
    ftsLst.propagation = fea.getPropagationCentroidQuad(voxi,voxr,muPix,ii,ftsLst.propagation,northDi,opts);
    %ftsLst.propagation = fea.getPropagationPixTracing(voxi,voxr,muPix,ii,ftsLst.propagation,northDi,0);
    %ftsLst.propagation = fea.getPropagationPixTracing(voxi,voxr,muPix,ii,ftsLst.propagation,northDi,1);
end

ftsLst.notes.propDirectionOrder = {'Anterior', 'Posterior', 'Left', 'Right'};

end




