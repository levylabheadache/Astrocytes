function [dFF, Fsub, baseF, C, A] = NormalizeFluor( Fastro, Fback, varargin )
%Convert raw fluorescence to dF/F, with background subtraction

IP = inputParser;
addRequired( IP, 'Fastro', @isstruct )
addRequired( IP, 'Fback', @iscell )
addParameter( IP, 'show', false, @islogical )
parse( IP, Fastro, Fback, varargin{:} ); 
show = IP.Results.show;

Nmovie = numel(Fback);
Ncell = size(Fastro(1).cell,2);
Nframe = cellfun(@numel, Fback);
% For each movie/astro/ROI subtract background region, then use OASIS to estimate base fluor, deconvolve Ca, activity. Then calculate dF/Fbase.
defaultStruct = repmat( struct('cell',[], 'soma',[], 'proc',[] ), 1, Nmovie );
Fsub = defaultStruct; C = defaultStruct; A = defaultStruct; baseF = defaultStruct; dFF = defaultStruct;
for m = flip(1:Nmovie)
    Fsub(m).cell = Fastro(m).cell - Fback{m}; 
    Fsub(m).soma = Fastro(m).soma - Fback{m}; 
    Fsub(m).proc = Fastro(m).proc - Fback{m}; 
    for c = flip(1:Ncell)
        [Ccell, Acell, Ocell] = deconvolveCa( Fsub(m).cell(:,c), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % Deconvolve dF/F  
        C(m).cell(:,c) = Ccell; 
        A(m).cell(:,c) = Acell; 
        baseF(m).cell(c) = Ocell.b; 
        [Csoma, Asoma, Osoma] = deconvolveCa( Fsub(m).soma(:,c), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); 
        C(m).soma(:,c) = Csoma; 
        A(m).soma(:,c) = Asoma; %
        baseF(m).soma(c) = Osoma.b; %
        [Cproc, Aproc, Oproc] = deconvolveCa( Fsub(m).proc(:,c), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); 
        C(m).proc(:,c) = Cproc; 
        A(m).proc(:,c) = Aproc; %
        baseF(m).proc(c) = Oproc.b; %
    end
    dFF(m).cell = (Fsub(m).cell - baseF(m).cell)./baseF(m).cell;
    dFF(m).soma = (Fsub(m).soma - baseF(m).soma)./baseF(m).soma;
    dFF(m).proc = (Fsub(m).proc - baseF(m).proc)./baseF(m).proc;
    if show
        close all;
        figure('WindowState','max', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        opt = {[0.05,0.08], [0.09,0.04], [0.06,0.06]};  % {[vert, horz], [bottom, top], [left, right] }
        for s = 5:-1:1, sp(s) = subtightplot(5,1,s, opt{:}); end
        subplot(sp(1));
        plot( Fastro(m).cell ); hold on;
        plot( Fback{m}, 'k' ); 
        title('Raw Fluorescence');
        subplot(sp(2));
        plot( Fsub(m).cell ); hold on;
        for c = 1:Ncell, line([1,Nframe(m)], baseF(m).cell(1,c)*[1,1], 'lineStyle','--', 'color','k' ); end
        title( sprintf('Backgroud Subtraction and Baseline (OASIS)' ) );
        subplot(sp(3));
        plot( dFF(m).cell );
        title('Normalization'); ylabel('dF/Fbase');
        subplot(sp(4));
        plot( C(m).cell );
        title('OASIS'); ylabel( '[Ca2+] (rel)' );
        subplot(sp(5));
        plot( A(m).cell );
        title('OASIS Deconvolution'); ylabel('Deconvolved Activity'); 
        linkaxes(sp,'x');
        xlim([1,Nframe(m)]);
        pause;
    end
end
end