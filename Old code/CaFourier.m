function [fSpec, Pspec, LFfrac, fGram, Tgram, Pgram] = CaFourier( T, fluor, varargin )
%CaFourier - Use Chronux toolbox to do fourier spectral analysis and calculate low frequency power

%if ispc, slm = '\'; else slm = '/'; end 
IP = inputParser;

addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'fluor', @isnumeric )
addParameter( IP, 'timeRes', 10, @isnumeric )
addParameter( IP, 'freqRes', 0.1, @isnumeric )
addParameter( IP, 'winInt', 3, @isnumeric )
addParameter( IP, 'lowCut', 0.2, @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, T, fluor, varargin{:} ); % mouse, comp, 
timeRes = IP.Results.timeRes; %timeRes = 20; % seconds 
freqRes = IP.Results.freqRes; %freqRes = 0.1; % Hz 
winInt = IP.Results.winInt; %winInt = 3; % seconds 
lowCut = IP.Results.lowCut;
show = IP.Results.show; 

Nroi = size(fluor,2);
TW = timeRes*freqRes;
K = 2*TW-1; 
chronuxParam.tapers = [TW, K]; 
chronuxParam.pad = 0;
chronuxParam.err = 0; % [2, 0.01];
chronuxParam.Fs = numel(T)/abs(T(end)-T(1)); %disp( chronuxParam );
movingWin = [timeRes, winInt ]; % [window length, interval between window centers] for spectrograms (optional)

fluor = fluor - mean(fluor,1); % Subtract the mean from each fluor channel
[Pspec, fSpec] = mtspectrumc( fluor, chronuxParam ); % Fourier power spectrum
LFfrac = sum(Pspec(fSpec <= lowCut,:),1)./sum( Pspec, 1 ); % low-frequency fraction
[Pgram,Tgram,fGram] = mtspecgramc( fluor,  movingWin, chronuxParam ); % spectrogram

if show
    opt = {[0.12,0.05], [0.1,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
    figure( 'WindowState','maximized','Color','w');
    for c = 1:Nroi
        % Fluorescence time series
        subtightplot(3,1,1, opt{:});  cla;
        plot( T, fluor(:,c) ); 
        title(sprintf('ROI %i', c));
        axis tight;
        % Spectrum
        subtightplot(3,1,2, opt{:});  cla;
        plot( fSpec, Pspec(:,c), 'LineWidth',3 ); hold on;
        xlabel('Frequency (Hz)'); ylabel('Power');
        title( sprintf('%2.3f pct power below %1.2f Hz', 100*LFfrac(c), lowCut) );
        set(gca,'Yscale','log');  % ,'FontSize',20
        axis tight;
        % Spectrogram
        subtightplot(3,1,3, opt{:});  cla;
        plot_matrix( Pgram(:,:,c), Tgram, fGram ); 
        xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('');
        colorbar('off')
        pause;
    end
end
