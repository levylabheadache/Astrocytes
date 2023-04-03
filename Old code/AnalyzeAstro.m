clear; clc; close all;
dataDir = 'D:\2photon\AB17\190530\FOV02\';

%[loadMovie, metadata, T, ~, speed] = LoadProcessed(dataDir, 'den');

MakeAstroROI( dataDir, '', true);
ReviewAstroROI( dataDir );