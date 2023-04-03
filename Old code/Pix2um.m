function p2uZoom = Pix2um(metadata)
%Pix2um calculates pixels per micron given magnification and objective

p2u16X = 0.53; %pixels per micron for our 16X objective, with 1X digital zoom 
p2uZoom = p2u16X*metadata.digiZoom; % pixels per micron for our 16X objective, with digital zoom 



end

