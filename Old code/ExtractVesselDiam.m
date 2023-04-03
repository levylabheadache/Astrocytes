function [bvDiam, lineDiam] = ExtractVesselDiam( vesselMovie, vessROI, varargin )
%
IP = inputParser;
addRequired( IP, 'vesselMovie', @isnumeric) 
addRequired( IP, 'vessROI', @isstruct)
addParameter( IP, 'Nthresh', 1, @isnumeric) 
addParameter( IP, 'show', false, @islogical )
parse( IP, vesselMovie, vessROI, varargin{:} ); 
show = IP.Results.show;

modelFcn = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) - b(4)./(1+exp(-(x-b(6))/b(5))) + b(7) ); % beta = [factor 1, growth rate 1, threshold 1, factor 2, growth rate 2, threshold 2, constant response offset]
Nvess = numel(vessROI);
maxDiam = nan(1,Nvess);
Nframe = size(vesselMovie,3);
gaussFilt = MakeGaussFilt( 3, 0, 1, 1, false ); % 3 frames wide, std dev = 1 pixel
lineDiam = cell(1, Nvess); bvDiam = nan(Nframe,Nvess); 
%figure('WindowState','max');
tic;
for b = 1:Nvess
    fprintf( sprintf('\n\n\n\nb = %i     ', b) );
    cropMovie = double( vesselMovie( vessROI(b).crop(2):vessROI(b).crop(2)+vessROI(b).crop(4), vessROI(b).crop(1):vessROI(b).crop(1)+vessROI(b).crop(3), : ) );
    maxDiam(b) = max( cellfun( @max, vessROI(b).lineDist ) );
    Nline = numel( vessROI(b).linePix ); % # of lines for this ROI
    tempR2 = nan(Nline, 1);  %tempResults = nan(Nline, 7); tempPval = nan(Nline, 7); 
    % Get the intensity profile along each line of the ROI and fit to a quasirectangular function
    [~, sortLineInd] = sort( vessROI(b).lineDiam, 'ascend' ); % sortLineDiam
    thinLineInd = sortLineInd(1:round(Nline/2))';
    lineValues = cell(1,Nline); smoothValues = cell(1,Nline); lineDiam{b} = nan(Nline, Nframe);
    for z = flip(1:Nframe)
        %fprintf( sprintf('[b,z] = [%i, %i]     ', b, z) );
        tempFrame = cropMovie(:,:,z);
        %{
        % Illustrate the process
        subplot(1,2,1); imshow( tempFrame, [] ); impixelinfo;
        subplot(1,2,2); 
        %}
        for p = thinLineInd %1:Nline
            lineValues{p}(:,z) = tempFrame( vessROI(b).lineInd{p} );
            smoothValues{p}(:,z) = filtfilt( gaussFilt, 1, lineValues{p}(:,z) ); 
            % Fit to step function model and read out the diameter
            modelFit = fitnlm( vessROI(b).lineDist{p}, smoothValues{p}(:,z), modelFcn, vessROI(b).fit.beta(p,:) );            
            %tempResults(p,:) = modelFit.Coefficients.Estimate;
            %tempPval(p,:) = modelFit.Coefficients.pValue;
            tempR2(p) =   modelFit.Rsquared.Adjusted;
            lineDiam{b}(p,z) = modelFit.Coefficients.Estimate(6) - modelFit.Coefficients.Estimate(3); %diff(vessROI(b).lineDist{p}( [maxPt, minPt] ));
            % Illustrate the process
            %{
            plot( vessROI(b).lineDist{p}, smoothValues{p}(:,z) ); hold on;
            xRange = linspace(0,vessROI(b).lineDist{p}(end))'; 
            plot( xRange, predict(modelFit, xRange), '--' );  % , 'color',[colorMat(p,:),0.5]
            title( sprintf('[b,z,p] = [%i, %i, %i]. Line Diam ~ %2.2f pix (R^2 = %2.2f)', b, z, p, lineDiam{b}(p,z), tempR2(p)  ) );
            pause;
            hold off;
            %}
        end
        estDiam = lineDiam{b}(:,z); 
        badEst = find( estDiam > maxDiam(b) | estDiam < 1 | tempR2 < 0.7 | isinf(tempR2) ); % suppress impossible/badly fit estimates of diameter
        %estDiam( estDiam > maxDiam(b) ) = NaN; estDiam( estDiam < 0 ) = NaN; estDiam( tempR2 < 0.8 ) = NaN; % suppress impossible/badly fit estimates of diameter
        if numel(badEst) < Nline
            estDiam(badEst) = NaN;
            if all(tempR2==0), error('all weights 0'); end
            bvDiam(z,b) = wmean( estDiam(estDiam>0),  tempR2(estDiam>0) );
        else
            warning('[b,z] = [%i, %i]  No good estimate of diameter obtained!', b, z);
        end
    end
    toc
end

if show
    clearvars sp;
    figure('WindowState','max');
    for b = flip(1:Nvess)
        sp(b) = subplot(Nvess, 1, b );
        plot( 1:Nframe, bvDiam(:,b) ); ylabel('Mean Diameter');
    end
    xlabel('Frame'); 
    linkaxes(sp,'x');
    xlim([-Inf,Inf]);
end
end

%{
    for n = 1:Nline
        
        
        %{
        tempMovie = false(  vessROI(b).cropSize(1), vessROI(b).cropSize(2), Nframe );
        tempMovie(vessROI(b).linePix{n}(:,2), vessROI(b).linePix{n}(:,1),: ) = true;
        XYZ = nan(Npix,2,Nframe);
        XYZ(:,1:2,:) = repmat( vessROI(b).linePix{n}, [1,1, Nframe] ); 
        XYZ(:,3,:) = repmat( 1:Nframe, [Npix, 1, 1] );
        
       pixInd = sub2ind(  [imSize(1), imSize(2), Nframe], XYZ )
        %}
        
    end
    %{
    % Crop the movie around the drawn ROI
    Xmin = min(vessROI(b).vertex(:,1)); Ymin = min(vessROI(b).vertex(:,2));
    cropWidth = max(vessROI(b).vertex(:,1)) - Xmin;
    cropHeight = max(vessROI(b).vertex(:,2)) - Ymin;
    vesselMovieCrop = vesselMean( Ymin:Ymin+cropHeight, Xmin:Xmin+cropWidth, : );
    vesselMeanCrop = imcrop( vesselMean, [Xmin, Ymin, cropWidth, cropHeight] ); 

    [cropSize(1), cropSize(2), ~] = size( vesselMovieCrop );
    %{
    vesselStdCrop = imcrop( vesselStd, [Xmin, Ymin, cropWidth, cropHeight] ); 
    vesselMinCrop = imcrop( vesselMin, [Xmin, Ymin, cropWidth, cropHeight] ); 
    vesselMaxCrop = imcrop( vesselMax, [Xmin, Ymin, cropWidth, cropHeight] ); 
    figure('WindowState','max'); 
    subplot(1,4,1); imshow( vesselMeanCrop, [] ); axis image;  title('Mean');
    subplot(1,4,2); imshow( vesselStdCrop, [] ); axis image; title('Std Dev');
    subplot(1,4,3); imshow( vesselMinCrop, [] ); axis image;  title('Min');
    subplot(1,4,4); imshow( vesselMaxCrop, [] ); axis image;  title('Max');
    impixelinfo; 
    pause;
    %}
    % Determine the lines that bisect the BV
    [sidePix{b,2}, ~] = lineSegPix( vessROI(b).vertex(2,:), vessROI(b).vertex(3,:), 0.4 ); 
    [sidePix{b,1}, ~] = lineSegPix( vessROI(b).vertex(1,:), vessROI(b).vertex(4,:), 0.4 );
    NsidePix(b,2) = size( sidePix{b,2}, 1 ); NsidePix(b,1) = size( sidePix{b,1}, 1 );
    cropPix{b,1} = sidePix{b,1} - repmat( [Xmin, Ymin], NsidePix(b,1), 1 ) + 1;
    cropPix{b,2} = sidePix{b,2} - repmat( [Xmin, Ymin], NsidePix(b,2), 1 ) + 1;

    meanThresh = multithresh( vesselMean(vessROI(b).ind), Nthresh);
    meanSeg = imquantize( vesselMeanCrop, meanThresh );
    meanSeg = meanSeg - 1;

    figure('WindowState','max'); 
    subplot(1,2,1); imshow( vesselMeanCrop, [] ); axis image;  title('Mean');
    subplot(1,2,2); imshow( meanSeg, [] ); axis image; title('Segmented'); hold on;
    plot( cropPix{b,1}(:,1), cropPix{b,1}(:,2), 'w.' );
    plot( cropPix{b,2}(:,1), cropPix{b,2}(:,2), 'w.' )
    cropPixSep = pdist2( cropPix{b,1}, cropPix{b,2} );
    goodLineInd = cell(NsidePix(b,1),1); goodLinePix = cell(NsidePix(b,1),1); goodLineDist = cell(NsidePix(b,1),1);
    for p = 1:NsidePix(b,1)
        linePix = cell(NsidePix(b,2),1); lineInd = cell(NsidePix(b,2),1); lineInt = cell(NsidePix(b,2),1); lineDist = cell(NsidePix(b,2),1); lineVessDiam = nan(NsidePix(b,2),1); 
        for n = 1:NsidePix(b,2)
            % Get pixel intensity and distance along the line across the BV
            linePix{n} = lineSegPix( cropPix{b,1}(p,:), cropPix{b,2}(n,:), 0.4 );
            lineInd{n} = sub2ind( cropSize, linePix{n}(:,2), linePix{n}(:,1) ); 
            lineDist{n} =  pdist2(cropPix{b,1}(p,:), linePix{n} )';
            lineInt{n} = meanSeg( lineInd{n} );
            lineVessInd = find( lineInt{n} );
            if numel(lineVessInd) > 1
                lineVessEdgeInd = lineVessInd([1,end]);
                lineVessDiam(n) = diff( lineDist{n}(lineVessEdgeInd) );
            end
            %plot( [cropPix{b,1}(p,1); cropPix{b,2}(n,1)], [cropPix{b,1}(p,2); cropPix{b,2}(n,2)], 'color',[0,0,0,0.5]);
            %title( sprintf('Line Diameter = %2.2f', lineVessDiam(n)) );
            %pause;
        end
        minDiam = min( lineVessDiam );
        nMin = find( lineVessDiam == minDiam, 1 ); 
        nGood = [];
        if numel(nMin) > 1 % if multiple lines have the same BV thickness, use the shortest pixel-to-pixel as a tie-breaker
            tempPixSep = cropPixSep(p,nMin);
            [~,tempInd] = min( tempPixSep );
            nGood = nMin( tempInd );
        elseif numel(nMin) == 1
            nGood = nMin;
        end
        if ~isempty( nGood )
            goodLineInd{p} = lineInd{nGood};
            goodLinePix{p} = linePix{nGood};
            goodLineDist{p} = lineDist{nGood};
            plot( [cropPix{b,1}(p,1); cropPix{b,2}(nGood,1)], [cropPix{b,1}(p,2); cropPix{b,2}(nGood,2)], 'color',[0,0,0,0.5]);
            title( sprintf('Line Diameter = %2.2f', lineVessDiam(nGood)) );
            pause;
        end
    end
    %}
    % What is the brightness along the side pixels?
    %{
    [sidePix{b,1}, ~] = lineSegPix( vessROI(b).vertex(1,:), vessROI(b).vertex(4,:), 0.4 );
    [sidePix{b,2}, ~] = lineSegPix( vessROI(b).vertex(2,:), vessROI(b).vertex(3,:), 0.4 );
    NsidePix(b) = size(sidePix{b,1}, 1);
    Nline(b) = size(sidePix{b,2}, 1);
    sideInd{b,1} = sub2ind( imSize, sidePix{b,1}(:,2), sidePix{b,1}(:,1) ); 
    sideInd{b,2} = sub2ind( imSize, sidePix{b,2}(:,2), sidePix{b,2}(:,1) ); 
    allSideInd = vertcat( sideInd{b,1}, sideInd{b,2} );
    Nind = numel(allSideInd);
    allSideInt = nan(Nind*Nframe,1);

    for z = 1:Nframe
         currentFrame = redMovie(:,:,z);
         allSideInt( (z-1)*Nind + (1:Nind) ) = currentFrame( allSideInd );
    end

    %redMaxThresh = 
    
    %colorMat = distinguishable_colors(NsidePix(b)); 
    edgeInd{b} = nan(NsidePix(b),2);  edgePix{b} = nan(NsidePix(b),2,2); 
    %edgeLine{b} =
    edgeDiam{b} = nan(NsidePix(b),1); % edgeAUC{b} = nan(NsidePix(b),1); 
    
     % For each pixel on one side of the vessel, find the nearest few pixels on the other side, draw line segments across, get the intensity profile under that line, estimate the BV width, and average 
    for p = 1:NsidePix(b)
        [~, tempSortInd] = sort( sidePixSep(p,:), 'ascend' );
        if show
            subplot(sp(1)); plot( sidePix{b,1}(p,1), sidePix{b,1}(p,2), 'wo', 'MarkerSize',6 );
        end
        linePix = cell(1,Nline); lineDist = cell(1,Nline); lineInt = cell(1,Nline);  
        bvEdgeInd = nan(Nline,2); lineDiam = nan(Nline,1); lineAUC = nan(Nline,1); 
        for n = 1:Nline 
            % Show the line
            %{
            subplot(sp(1))
            plot( sidePix{b,2}(tempSortInd(n),1), sidePix{b,2}(tempSortInd(n),2), 'x', 'MarkerSize',10 ); % , 'color',colorMat(n,:)
            plot( [sidePix{b,1}(p,1); sidePix{b,2}(tempSortInd(n),1)], [sidePix{b,1}(p,2); sidePix{b,2}(tempSortInd(n),2)], 'MarkerSize',10, 'Color',[colorMat(n,:),alphaValue] );
            %}
            % Get pixel intensity and distance along the line across the BV
            linePix{n} = lineSegPix( sidePix{b,1}(p,:), sidePix{b,2}(tempSortInd(n),:), 0.4 );
            lineInd = sub2ind( imSize, linePix{n}(:,2), linePix{n}(:,1) ); 
            lineDist{n} =  pdist2(sidePix{b,1}(p,:), linePix{n} )';
            lineInt{n} = currentFrame( lineInd );
            bvInd = find(lineInt{n} >= bvThresh(z,b));
            if ~isempty(bvInd)
                bvEdgeInd(n,:) = lineInd( bvInd([1,end]) )';
                lineDiam(n) = diff( lineDist{n}(bvInd([1,end])) ); % lineDist{n}(bvInd(end))-lineDist{n}(bvInd(1));
                %lineAUC(n) = sum(lineInt{n});
            end
            % Plot the current line's intensity profile
            %{
            plot( lineDist{n},  lineInt{n}, 'Color', [0,0,0,0.5] ); hold all;
            line( [lineDist{n}(1), lineDist{n}(end)], bvThresh(z,b)*[1,1], 'lineStyle','--','color',[0,0,0,0.5] )
            xlabel('Distance Along the Line'); ylabel('Pixel Intensity'); title( sprintf('diameter = %2.2f', lineDiam(n) ) );
            %}
        end
    
    
    % {
    Xmin = min(vessROI(b).vertex(:,1)); Ymin = min(vessROI(b).vertex(:,2));
    cropWidth = max(vessROI(b).vertex(:,1)) - Xmin;
    cropHeight = max(vessROI(b).vertex(:,2)) - Ymin;
    [sidePix{b,1}, ~] = lineSegPix( vessROI(b).vertex(1,:), vessROI(b).vertex(4,:), 0.4 );
    [sidePix{b,2}, ~] = lineSegPix( vessROI(b).vertex(2,:), vessROI(b).vertex(3,:), 0.4 );
    sidePixSep = pdist2( sidePix{b,1}, sidePix{b,2} );
    NsidePix(b)= size(sidePix{b,1}, 1);
    colorMat = distinguishable_colors(NsidePix(b)); 
    edgeInd{b} = nan(NsidePix(b),2);  edgePix{b} = nan(NsidePix(b),2,2); edgeDiam{b} = nan(NsidePix(b),1); % edgeAUC{b} = nan(NsidePix(b),1); 
    for z = 1:Nframe
        currentFrame = redMovie(:,:,z);
        if show
            subplot(sp(1)); imshow( currentFrame, [], 'InitialMagnification','fit' ); hold all; axis off; impixelinfo;
            subplot(sp(2)); cla; 
        end
        % Determine intensity threshold
        cropIm = imcrop( currentFrame, [Xmin, Ymin, cropWidth, cropHeight] ); 
        bvThresh(z,b) = double(multithresh( cropIm ));
        %{
        figure;
        imshow( cropIm, [] );
        title( sprintf('threshold = %i', bvThresh(z,b) ) );
        impixelinfo;
        pause;
        %}
        % For each pixel on one side of the vessel, find the nearest few pixels on the other side, draw line segments across, get the intensity profile under that line, estimate the BV width, and average 
        for p = 1:NsidePix(b)
            [~, tempSortInd] = sort( sidePixSep(p,:), 'ascend' );
            if show
                subplot(sp(1)); plot( sidePix{b,1}(p,1), sidePix{b,1}(p,2), 'wo', 'MarkerSize',6 );
            end
            linePix = cell(1,Nline); lineDist = cell(1,Nline); lineInt = cell(1,Nline);  
            bvEdgeInd = nan(Nline,2); lineDiam = nan(Nline,1); lineAUC = nan(Nline,1); 
            for n = 1:Nline 
                % Show the line
                %{
                subplot(sp(1))
                plot( sidePix{b,2}(tempSortInd(n),1), sidePix{b,2}(tempSortInd(n),2), 'x', 'MarkerSize',10 ); % , 'color',colorMat(n,:)
                plot( [sidePix{b,1}(p,1); sidePix{b,2}(tempSortInd(n),1)], [sidePix{b,1}(p,2); sidePix{b,2}(tempSortInd(n),2)], 'MarkerSize',10, 'Color',[colorMat(n,:),alphaValue] );
                %}
                % Get pixel intensity and distance along the line across the BV
                linePix{n} = lineSegPix( sidePix{b,1}(p,:), sidePix{b,2}(tempSortInd(n),:), 0.4 );
                lineInd = sub2ind( imSize, linePix{n}(:,2), linePix{n}(:,1) ); 
                lineDist{n} =  pdist2(sidePix{b,1}(p,:), linePix{n} )';
                lineInt{n} = currentFrame( lineInd );
                bvInd = find(lineInt{n} >= bvThresh(z,b));
                if ~isempty(bvInd)
                    bvEdgeInd(n,:) = lineInd( bvInd([1,end]) )';
                    lineDiam(n) = diff( lineDist{n}(bvInd([1,end])) ); % lineDist{n}(bvInd(end))-lineDist{n}(bvInd(1));
                    %lineAUC(n) = sum(lineInt{n});
                end
                % Plot the current line's intensity profile
                %{
                plot( lineDist{n},  lineInt{n}, 'Color', [0,0,0,0.5] ); hold all;
                line( [lineDist{n}(1), lineDist{n}(end)], bvThresh(z,b)*[1,1], 'lineStyle','--','color',[0,0,0,0.5] )
                xlabel('Distance Along the Line'); ylabel('Pixel Intensity'); title( sprintf('diameter = %2.2f', lineDiam(n) ) );
                %}
            end
            [~, nMin] = min( lineDiam );
            edgeInd{b}(p,:) = bvEdgeInd(nMin,:);
            [Ytemp, Xtemp] = ind2sub( imSize, edgeInd{b}(p,:) );
            edgePix{b}(p,:,2) = [Xtemp(2), Ytemp(2)];  edgePix{b}(p,:,1) = [Xtemp(1), Ytemp(1)];
            edgeDiam{b}(p) = lineDiam(nMin);
            %edgeAUC{b}(p) = lineAUC(nMin);
            bvDiam(z,b) = mean( edgeDiam{b} );
            if show
                % Plot the shortest line's intensity profile
                subplot(sp(2));
                plot( lineDist{nMin},  lineInt{nMin}, 'Color', [colorMat(p,:),0.5] ); hold all;
                line( [lineDist{nMin}(1), lineDist{nMin}(end)], bvThresh(z,b)*[1,1], 'lineStyle','--','color',[0,0,0,0.5] )
                % Plot the shortest line
                subplot(sp(1));
                plot( sidePix{b,2}(tempSortInd(nMin),1), sidePix{b,2}(tempSortInd(nMin),2), 'x', 'MarkerSize',10 ); hold on;
                plot( [sidePix{b,1}(p,1); sidePix{b,2}(tempSortInd(nMin),1)], [sidePix{b,1}(p,2); sidePix{b,2}(tempSortInd(nMin),2)], 'MarkerSize',10, 'Color',[colorMat(p,:),alphaValue] );
                plot( [edgePix{b}(p,1,1); edgePix{b}(p,1,2)], [edgePix{b}(p,2,1); edgePix{b}(p,2,2)], '.','MarkerSize',10, 'Color',colorMat(p,:) )
            end
        end
        if show
            title( sprintf('[z,b] = [%i, %i]  diameter = %2.2f', z,b, bvDiam(z,b) ) );
            pause(0.2);
        end
    end
    toc;
    %}
%}
   %{
    smoothDeriv = cellfun( @diff, smoothValues, 'UniformOutput',false );
    [~, minInds] = cellfun( @min, smoothDeriv, 'UniformOutput',false  );
    [~, maxInds] = cellfun( @max, smoothDeriv, 'UniformOutput',false  );
    
    lineDiam{b} = nan(Nline, Nframe);
    %figure;
    for z = 1:Nframe
        for p = 1:Nline
            minPt = minInds{p}(z)+1;  maxPt = maxInds{p}(z);
            lineDiam{b}(p,z) =  diff(vessROI(b).lineDist{p}( [maxPt, minPt] ));
            %{
            subplot(2,1,1); cla;
            plot( vessROI(b).lineDist{n}, smoothValues{n}(:,z) ); hold on;
            plot( vessROI(b).lineDist{n}( minPt), smoothValues{n}(minPt,z), 'o' );
            plot( vessROI(b).lineDist{n}( maxPt), smoothValues{n}(maxPt,z), 'o' );
            ylabel('Intensity');
            title( sprintf('[b,z,n] = [%i, %i, %i]. Estimated Diameter = %2.2f pixels', b, z, n, lineDiam(n,z) ) );

            subplot(2,1,2); cla;
            plot( vessROI(b).lineDist{n}(2:end), smoothDeriv{n}(:,z) ); hold on;
            plot( vessROI(b).lineDist{n}( minPt), smoothDeriv{n}(minInds{n}(z),z), 'o' );
            plot( vessROI(b).lineDist{n}( maxPt+1), smoothDeriv{n}(maxInds{n}(z),z), 'o' );
            ylabel('\Delta Intensity');
            pause;
            %}
        end
       
    end
    lineDiam{b}( lineDiam{b} <= 1 ) = NaN; % imagesc( lineDiam{b} ); colorbar;
    bvDiam(:,b) = mean( lineDiam{b}, 1, 'omitnan' )';
    %}
    %{
    if show
        colorMat = distinguishable_colors(Nline);
        figure; sp(2) = subplot(1,2,2); sp(1) = subplot(1,2,1); 
        for z = flip(1:Nframe)
            subplot(sp(2)); cla;
            subplot(sp(1)); cla;
            imshow( tempFrame, [] ); hold on;
            title( sprintf('b = %i, z = %i', b, z) );
            for n = 1:Nline
                subplot(sp(1))
                plot( vessROI(b).linePix{n}(:,1), vessROI(b).linePix{n}(:,2), 'x', 'color', colorMat(n,:) )
                subplot(sp(2))
                plot( vessROI(b).lineDist{n}, smoothValues{n}(:,z), 'color',[colorMat(n,:),0.5] ); hold on;
            end
            pause(0.1);
        end
    end
    %}