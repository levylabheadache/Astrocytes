function LFratio = LowFreqRatio(fluorData, chronuxParam, lowFreq) 

fluorData = fluorData - mean(fluorData,1); % subtract DC offset
[Stemp, tempFreq] = mtspectrumc( fluorData, chronuxParam ); % get Fourier power spectrum

%figure;
%subplot(2,1,1); plot( fluorData );
%subplot(2,1,2); plot( tempFreq, Stemp );

LFframes = find(tempFreq > lowFreq(1) & tempFreq <= lowFreq(2));
LFratio = abs( sum( Stemp(LFframes,:), 1 )./sum(Stemp, 1) );  %#ok<FNDSB> % ( find(tempFreq >= 0 & tempFreq <= 1)


end

