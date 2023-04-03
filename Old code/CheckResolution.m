clear
mouse = 'AB23';
date = '191017';
fov = 1;


[~, dataDir] = SetDataPath(mouse, date, [], fov);

metadata = LoadProcessed( dataDir, 'pool' );
spatialRes = 2/Pix2um( metadata(1) ) %factor of 2 due to pixel-binning for denoising

temporalRes = 3/15.49 % factor of 3 due to frame-binning for denoising