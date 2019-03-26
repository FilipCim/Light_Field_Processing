clear all
close all
clc
%%
utilities_functions = './utilities';
mps_functions = './utilities/mps';
addpath(utilities_functions, mps_functions); 

% SeDuMi path, used for initial testing (already precompiled)
addpath('C:\development\SeDuMi_1_3\');
% SPAMS path, used for final reconstruction
SPAMSPath = 'C:\development\spams-matlab-v2.6_Copy';
addpath(SPAMSPath);

currentDirectory = 'D:\faks\projekt\script\GitRepository';

closescreen

%% ACQUISITION PARAMETERS
% camera and projector settings for slideshow acquisition

% delay time in seconds
params.delayTimeSlideShow = 0.5;
params.delayTimeCamera = 0.5;
params.blockSize = 8;
params.projectorResolutionDecimationFactor = 2;

%%%%%%%%%%%%% projector settings %%%%%%%%%%%%%%
% screenID where you want to project
params.screenID = 2;
% size of the screen
params.screenWidth = 1920;
params.screenHeight = 1080;
% screen bitdepth
params.bitDepth = 8;

% %%%%%%%%%%%%% camera settings %%%%%%%%%%%%%
% cam = gigecam;
% cam.ExposureTime = 8333.333333333;
% cam.AcquisitionFrameRate = 10;
% cam.AcquisitionFrameRateAbs = 10;
% cam.GainAutoBalance = 'off';
% cam.BlackLevelAutoBalance = 'off';
% cam.AcquisitionBurstFrameCount = 1;
% cam.GevSCPSPacketSize = 4000;
% params.roi(2) = 0;
% params.roi(4) = 2058;
%
% params.cam = cam;
params.writeFlag = 1;

debug = 1;
% %%
% if debug
% % Loading unpacked Light Field (here so the unpacking wouldn't be done
% % every time ROI is calculated)
% ROI_image = load('.\data\Scenes\ROI_image2Kosa.mat');
% img = ROI_image.img;
% roi = ROI_f(img, debug);
% 
% % Increasing ROI for 20% and checking the size of the ROI so it would not
% % extend the size of the image
% 
% % 10% of the width (x dim)
% roi_scale(1) = roi(3)*0.1;
% % 10% of the length (y dim)
% roi_scale(2) = roi(4)*0.1;
% 
% roi_wide(1) = round(roi(1) - roi_scale(1));
% if roi_wide(1) < 0
%     roi_wide(1) = 1;
% end
% 
% roi_wide(2) = round(roi(2) - roi_scale(1));
% if roi_wide(2) < 0
%     roi_wide(2) = 1;
% end
% 
% roi_wide(3) = round(roi(3) + 2*roi_scale(1));
% if (roi_wide(1)+roi_wide(3)) > size(img, 2)
%     roi_wide(3) = size(img, 1) - roi_wide(1);
% end
% 
% roi_wide(4) = round(roi(4) + 2*roi_scale(2));
% if (roi_wide(4)+roi_wide(2)) > size(img, 1)
%     roi_wide(4) = size(img, 1) - roi_wide(2);
% end
% 
% if debug
%         figure(1122)
%         imshow(img);
%         hold on;
%         rectangle('Position',...
%             [roi_wide(1),roi_wide(2),roi_wide(3),roi_wide(4)],...
%             'EdgeColor','r','LineWidth',2 )
%         drawnow
%         hold off;
%     end
% 
% disp('Found ROI...')
% end
%% CALCULATE PER PIXEL CORRESPONDENCES USING MPS FOR CALIBRATION

imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\mpsImagesTest\');
mpsDirectory = imagesCameraDirectoryPath;

[rows, cols, cameraPoints, projectorPoints, ~] = decodeMPS(mpsDirectory, params);

% rowsInt = round(rows);
% colsInt = round(cols);

[~,~,rowsBlock] = histcounts(round(rows), round((params.screenHeight/params.blockSize)/params.projectorResolutionDecimationFactor));
[~,~,colsBlock] = histcounts(round(cols), round((params.screenWidth/params.blockSize)/params.projectorResolutionDecimationFactor));


%% CALCULATE HOMOGRAPHY BETWEEN CAMERA AND PROJECTOR


% tformCamera2Projector = fitgeotrans(cameraPoints, projectorPoints, 'projective');
% outputView = imref2d([params.screenHeight, params.screenWidth]);

tformProjector2Camera = fitgeotrans(projectorPoints, cameraPoints, 'projective');
outputView = imref2d([size(rows,1), size(rows,2)]);

% colsBlock = round(imwarp(colsBlock, tformCamera2Projector, 'OutputView',outputView));
% rowsBlock = round(imwarp(rowsBlock, tformCamera2Projector, 'OutputView',outputView));

%% ESTIMATE VIGNETTING FUNCTION

imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\whiteImagesTest\');

dirList = dir(imagesCameraDirectoryPath);
isFile = ~[dirList.isdir];
imageFilenames = sort({dirList(isFile).name});


for i = 1:length(imageFilenames)
    imagesCamera(:,:,i) = imread([imagesCameraDirectoryPath, imageFilenames{i}]);
%     figure(100)
%     imagesc(imagesCamera(:,:,i));
%     waitforbuttonpress
%     axis image
end

whiteImageCamera= mean(imagesCamera,3);
ProjectorImageWarped = imwarp(ones([params.screenHeight, params.screenWidth]), tformProjector2Camera, 'OutputView',outputView);

whiteImageCamera = double(whiteImageCamera)./double(max(whiteImageCamera(:)));

vignettingFunctionEstimation = double(ProjectorImageWarped)./whiteImageCamera;

figure,
imagesc(vignettingFunctionEstimation)
title('Estimated vignetting function')

%%
% Importing white image...
imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\mpsImagesTest\');

whiteImageScene = imread([imagesCameraDirectoryPath, 'white', '.tiff']);
whiteImageScene = double(whiteImageScene)./double(max(whiteImageScene(:)));

figure
imagesc(whiteImageScene.*vignettingFunctionEstimation)

%% ESTIMATE NOISE

imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\blackImagesTest\');

dirList = dir(imagesCameraDirectoryPath);
isFile = ~[dirList.isdir];
imageFilenames = sort({dirList(isFile).name});


for i = 1:length(imageFilenames)
    imagesCamera(:,:,i) = imread([imagesCameraDirectoryPath, imageFilenames{i}]);
    
    imagesCameraNoise(:,:,i) = imagesCamera(:,:,i);
    
end

noiseEstimation = mean(imagesCameraNoise, 3).*vignettingFunctionEstimation;

figure,
imagesc(noiseEstimation)
title('Estimated noise')


%% CALCULATE PER PIXEL CORRESPONDENCES USING MPS FOR SCENE

imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\FerScene\mpsImagesTest\');

[rows, cols, cameraPoints, projectorPoints, ~] = decodeMPS(imagesCameraDirectoryPath, params);

% rowsInt = round(rows);
% colsInt = round(cols);

[~,~,rowsBlock] = histcounts(round(rows), round((params.screenHeight/params.blockSize)/params.projectorResolutionDecimationFactor));
[~,~,colsBlock] = histcounts(round(cols), round((params.screenWidth/params.blockSize)/params.projectorResolutionDecimationFactor));

% rows - a matrix that shows which index of the projector corresponds to
% some pixel on camera (because we are working with a image from the
% camera. Indexes are sorted like it's a vector (reading the matrix top
% down, left right). This way we can see that projector pixel with the
% index R (we can get the index from Phase Shift coding scheme we used)
% corresponds to the X,Y pixel on the camera. Here we can see that Phase
% Shifting was used to be able to get the sub-pixel correspondencies (in
% between the maxs and mins of sinuses we interpolate sine function and get
% the exact index (with the decimal point) that we need

% cols - same as rows, but for columns (vertical projections of sine images
% on the scene)

% cameraPoints - pairs of points on camera and projector that have good
% correspondency. These points are not corresponding to 2 points or half of
% point or something similar. These can be used for homography.

% projectorPoints - same as cameraPoints but the other way around

% rowsBlock - we are using blocks in our processing so we need a value for
% each block in the image. By using histcounts we can get the values in
% bins for every block

% colsBlock - same as rowBlock, but for columns


%% SELECT ROI FOR CS RECONSTRUCTION 
firstPhaseImage = ((imread(fullfile(currentDirectory,... 
    '\data\Scenes\FerScene\hadamardImagesTest\measurement_mask_1_01.tiff'))));
%firstPhaseImage = double(firstPhaseImage);

figure(111),
imshowpair(whiteImageScene, firstPhaseImage)

[xi, yi] = getpts(gcf);
close gcf;

position(1)=xi(1);
position(2)=yi(1);
position(3)=xi(2)-xi(1);
position(4)=yi(2)-yi(1);

% ROI selection done but manually inputing two points on the image. First
% one is upper left. Second one is lower right. Last point needs to be
% double-click to be able to add the last point. 

% The same can be done with not casting to double and pressing enter after
% selecting the last point.

startRow = (impixel(rowsBlock, xi(1), yi(1)));
startRow = startRow(1);
endRow = impixel(rowsBlock, xi(2), yi(2));
endRow = endRow(1);
startCol = impixel(colsBlock, xi(1), yi(1));
startCol = startCol(1);
endCol = impixel(colsBlock, xi(2), yi(2));
endCol = endCol(1);

%%
% This is liner sorting of the indexes for block correspondencies.
blockCorrespondences = (rowsBlock-1).*max(colsBlock(:))+colsBlock;

% phasesMask{1} = (mod(rowsBlock,2) & mod(colsBlock,2));
% phasesMask{2} = (mod(rowsBlock,2) & mod(colsBlock,2)==0);
% phasesMask{3} = (mod(rowsBlock,2)==0 & mod(colsBlock,2));
% phasesMask{4} = (mod(rowsBlock,2)==0 & mod(colsBlock,2)==0);

phasesMask{1} = (mod(rowsBlock,2) & mod(colsBlock,2));
phasesMask{2} = (mod(rowsBlock,2) & mod(colsBlock,2)==0);
phasesMask{3} = (mod(rowsBlock,2)==0 & mod(colsBlock,2));
phasesMask{4} = (mod(rowsBlock,2)==0 & mod(colsBlock,2)==0);


for phaseNo = 1:4
    phaseNo
    
    % crop a phase from the image by using phasesMask
    phases{phaseNo} = (blockCorrespondences .* phasesMask{phaseNo});
   
    % using watershed agorithm to get a one pixel border between two
    % blocks in one phase
    phasesWatershed{phaseNo} = (watershed(-bwdist(~phases{phaseNo}))>0);
    
    s = regionprops(phasesMask{phaseNo}, 'Area', 'PixelList', 'PixelIdxList', 'Centroid');
    % median filtering - set 0 to all of the small values (20% of the
    % median of the area
    idx = cell2mat({s.Area}) < 0.2*median(cell2mat({s.Area}));
    s(idx) = [];
    
    sWater = regionprops(phasesWatershed{phaseNo}, 'Area', 'PixelList', 'PixelIdxList','Centroid');
    
    centroidCoords = cell2mat({s.Centroid}');
    centroidCoordsWater = cell2mat({sWater.Centroid}');
    
    % euclidean distance between centroids of the watershed and
    % non-watershed version of the phase
    D = pdist2(cell2mat({sWater.Centroid}'), cell2mat({s.Centroid}'), 'euclidean');
    
    [val, idx] = min(D,[],1);
    
    phasesWatershed{phaseNo} = double(phasesWatershed{phaseNo});
    
    % ????????????????????????????????????????????????????????????????????
    for i = 1:size(idx,2)
        pixValue = repmat(impixel(phases{phaseNo}, centroidCoords(i,1),centroidCoords(i,2)), size(phasesWatershed{phaseNo}(sWater(idx(i)).PixelIdxList)));
        
        phases{phaseNo}(sWater(idx(i)).PixelIdxList) = pixValue(:,1);
    end
    
    
end




%%

% show all of the markers and bordes for phases
figure,
subplot(221), imshowpair(phases{1}, firstPhaseImage), axis image
hold on
plot(xi(1), yi(1), 'rx', 'MarkerSize', 10)
plot(xi(2), yi(2), 'rx', 'MarkerSize', 10)
subplot(222), imagesc(phases{2}), axis image
hold on
plot(xi(1), yi(1), 'rx', 'MarkerSize', 10)
plot(xi(2), yi(2), 'rx', 'MarkerSize', 10)
subplot(223), imagesc(phases{3}), axis image
hold on
plot(xi(1), yi(1), 'rx', 'MarkerSize', 10)
plot(xi(2), yi(2), 'rx', 'MarkerSize', 10)
subplot(224), imagesc(phases{4}), axis image
hold on
plot(xi(1), yi(1), 'rx', 'MarkerSize', 10)
plot(xi(2), yi(2), 'rx', 'MarkerSize', 10)


%% Measuring
imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\FerScene\hadamardImagesTest\');

% processing phase-by-phase
% Here we iterate every phase for every measurement (in the hadamard 
% folder)
for phaseNumber = 1:4
    for measurementNumber = 1:64
        display(['Processing:     (',num2str((phaseNumber-1)*64+measurementNumber),'/',num2str(4*64),')', imagesCameraDirectoryPath, 'measurement_mask_', num2str(phaseNumber), '_', num2str(measurementNumber, '%02d'), '.tiff'])
        
        image = imread([imagesCameraDirectoryPath, 'measurement_mask_', num2str(phaseNumber), '_', num2str(measurementNumber, '%02d'), '.tiff']);
        image = double(image);
        
        % Noise reduction
        image = double(image).*vignettingFunctionEstimation;
        image = image - noiseEstimation;   
        
        % Setting all NaN values to 0, NaNs caused by devignetting,
        % dividing by pixels with value 0 in white image (edges)
        image(isnan(image)) = 0;
        
        tic
        
        % depending on what is the current phase that is being processed,
        % this will choose it 
        switch mod(phaseNumber,4)
            case 1
                
                % every two rows, every two colums
                for row = startRow:2:endRow
                    for col = startCol:2:endCol
                        
                        % take the indexes and sort them linearly
                        % (vectorize them)
                        
                        mask = (phases{1}==((row-1)*max(colsBlock(:))+col));
                        %mask = (phases{1}==((row+startRow-1)*(max(colsBlock(:))-startCol) +col+startCol));
                        
                        % take the right values and just calculate a mean
                        % value of the whole block. That's it.
                        measurements(row-startRow+1, col-startCol+1, measurementNumber)= mean2((image(mask)));
                        
%                         figure(101)
%                         imshowpair(image, mask)
%                         drawnow
                        
                    end
                end
                
            case 2
                for row = startRow:2:endRow
                    for col = startCol+1:2:endCol
                        
                        mask = (phases{2}==((row-1)*max(colsBlock(:))+col));
                        % mask = (phases{2}==((row+startRow-1)*max(colsBlock(:)) +col+startCol));
                        measurements(row-startRow+1, col-startCol+1, measurementNumber)= mean2((image(mask)));
                        
%                         figure(101)
%                         imshowpair(image, mask)
%                         drawnow

                        
                    end
                end
                
            case 3
                for row = startRow+1:2:endRow
                    for col = startCol:2:endCol
                        
                        mask = (phases{3}==((row-1)*max(colsBlock(:))+col));
                        % mask = (phases{3}==((row+startRow-1)*max(colsBlock(:)) +col+startCol));
                        measurements(row-startRow+1, col-startCol+1, measurementNumber)= mean2((image(mask)));

%                         figure(101)
%                         imshowpair(image, mask)
%                         drawnow
                        
                    end
                end
                
            case 0
                for row = startRow+1:2:endRow
                    for col = startCol+1:2:endCol
                        
                        mask = (phases{4}==((row-1)*max(colsBlock(:))+col));
                        % mask = (phases{4}==((row+startRow-1)*max(colsBlock(:)) +col+startCol));
                        measurements(row-startRow+1, col-startCol+1, measurementNumber)= mean2((image(mask)));
% 
%                         figure(101)
%                         imshowpair(image, mask)
%                         drawnow
                        
                    end
                end
                
   
        end
                
        toc
        
        
    end
end

display('Measurement preprocessing ended!')

%% Check measurments
figure()
colormap gray
imagesc(measurements(:,:,22));
title('One of the measurements');

%% CREATE MEASUREMENT MATRIX FROM PROJECTED MEASUREMENT MASKS

[synth_masks, synth_mask_number_of_ones, phi] = loadMeasurementMasks... 
    ('.\data\imagesProjector\hadamardMeasurementMaskImages\Mask Tiles_8x8\', '.png');

%% MEASUREMENT LINEARIZATION
% measurement linearization process - masks with different percentage of ones were
% projected on white and black paper (number of ones on white paper is n and
% number of ones on black paper is 33-n) - this corresponds to cases when
% different percentage of bright pixels is projected on white or on black
% part of measured scene
% sum of measurements obtained on white and black paper is calculated
% the sum is used to linearize measurements

% corrected_measurements = linearizeMeasurements(measurementsS, measurementsW, measurementsB);

block_size = 8;

corrected_measurements = linearizeMeasurements(measurements, measurements, measurements);
% corrected_measurements = corrected_measurements .* repmat(reshape(synth_mask_number_of_ones, [1,1,size(measurementsS,3)]), [size(measurementsS,1), size(measurementsS,2), 1]);
% net
% corrected_measurements = corrected_measurements + 0.1*rand(size(measurementsW,1), size(measurementsW,2), 64);

% for i = 1:64
% %     corrected_measurements(:,:,i) = measurementsS(:,:,i)./calibMeasurements1(:,:,i);
%     corrected_measurements(:,:,i) = measurementsS(:,:,i)./calibMeasurements2(:,:,i);
%
% end

% noOfMeasurementsReconstruction = ceil(measurementRate*block_size^2);

% y = corrected_measurements;
y = measurements;
ind = 2:30;

y=y(:,:,ind);
% reduce measurement matrix to defined number of measurements
phi_r = phi(ind, :)./sum(phi(ind, :), 2);


% define transformation matrix type(dct or dwt) and if dwt selected define
% wavelet type used to generate dwt transformation matrix
% psi_inv = dctmtx(block_size);
psi_inv = wmpdictionary(block_size ,'lstcpt',{'dct'});
% psi_inv = wmpdictionary(block_size ,'lstcpt',{'dct','RnIdent', 'sin', 'poly', {'sym4',5}, {'sym4',2}, {'sym4',3}, {'sym4',7},{'wpsym4',5}, {'wpsym4',10} , {'wpsym4',8}});
% psi_inv = wmpdictionary(block_size ,'lstcpt',{{'haar',3}});

psi_inv=kron(psi_inv, psi_inv);
% matrix theta definition: y = theta * x

theta = sparse(phi_r * psi_inv); % phi_r * psi^(-1)


% define reconstruction package parameters
reconstructionOptions.Package = 'mex';
% for cvx reconstruction package user has choice of algorithms(sedumi, sdpt3)
% reconstructionOptions.Algorithm = 'sedumi';
reconstructionOptions.lambda = 0;

reconstruction=[];

if strcmp(lower(reconstructionOptions.Package), 'mex')
    oldPath = cd(SPAMSPath);
    start_spams;
    cd(oldPath);
end


for i = 1:size(measurements, 1)
    for j = 1:size(measurements, 2)
        
        subimage_estimation = l1Optimization(y(i, j, :), psi_inv, theta, reconstructionOptions);
        
        subimage_estimation = reshape(subimage_estimation, block_size, block_size);
        reconstruction((i-1)*block_size+1:(i-1)*block_size+block_size,(j-1)*block_size+1:(j-1)*block_size+block_size) = subimage_estimation;
        
    end
    
    figure(101)
    imagesc(reconstruction), axis image, colormap gray
    drawnow
end

%% CALCULATE CALIBRATION PATCHES


coeffOfVar = std(measurements./repmat(mean(measurements,3), [1,1,64]),0,3);
coeffOfVar = coeffOfVar./max(coeffOfVar(:));

% figure(), imagesc(coeffOfVar)


% imagesc(coeffOfVar<0.975*median(coeffOfVar(:)))

%% Scaling reconstruction
reconstruction = mat2gray(reconstruction);

waitforbuttonpress

%% INVERSE HOMOGRAPHY
% All of measurements and reconstruction are made in the Projector plane.
% With homography (inverse one) plane of viewing is shifted to the Camera
% plane.
% This should be done for every PixelView image in the tiled format.
% Reconstruction is the same, only thing that is different is homography,
% so MPS should be calculated for every image, but reconstruction shouldn't
debug = 1;
homograph = 0;

% Prealocating cells 
reconstruction_cell = cell(15,15);
homography_cell = cell(15,15);
image_cell = cell(15,15);

% Variable homograph dictates whether or not the (inverse) homography is
% calculated
if homograph
    tic
% Calculating ROI 
% Loading unpacked Light Field (here so the unpacking wouldn't be done
% every time ROI is calculated)
ROI_image = load('.\data\Scenes\ROI_image2Kosa.mat');
img = ROI_image.img;
roi = ROI_f(img, debug);

% Increasing ROI for 20% and checking the size of the ROI so it would not
% extend the size of the image

% 10% of the width (x dim)
roi_scale(1) = roi(3)*0.1;
% 10% of the length (y dim)
roi_scale(2) = roi(4)*0.1;

roi_wide(1) = round(roi(1) - roi_scale(1));
if roi_wide(1) < 0
    roi_wide(1) = 1;
end

roi_wide(2) = round(roi(2) - roi_scale(1));
if roi_wide(2) < 0
    roi_wide(2) = 1;
end

roi_wide(3) = round(roi(3) + 2*roi_scale(1));
if (roi_wide(1)+roi_wide(3)) > size(img, 2)
    roi_wide(3) = size(img, 1) - roi_wide(1);
end

roi_wide(4) = round(roi(4) + 2*roi_scale(2));
if (roi_wide(4)+roi_wide(2)) > size(img, 1)
    roi_wide(4) = size(img, 1) - roi_wide(2);
end

if debug
        figure(112)
        imshow(img);
        hold on;
        rectangle('Position',...
            [roi_wide(1),roi_wide(2),roi_wide(3),roi_wide(4)],...
            'EdgeColor','r','LineWidth',2 )
        drawnow
        hold off;
end

disp('Found ROI...')

% Calculating homography for all images in tiled format with skipping the
% corners. Criteria for skipping is the sum of values of all pixels. The
% threshold is set experimentally (with inspecting images in the corners).
for i = 1:15
    for j = 1:15
        tic
        
        TiledImagePath = '.\data\Scenes\FerScene\mpsImagesTiled\';
        SaveDirPath = '.\data\Scenes\FerScene\mpsImagesTest\';

        dirList = dir(TiledImagePath);
        isFile = ~[dirList.isdir];
        imageFilenames = {dirList(isFile).name};

        % For all MPS images and for one pixel view (one tile in tiled image)
        for g = 1:size(imageFilenames, 2)
            image = imread([TiledImagePath, imageFilenames{1,g}]);
            cutImage = CutTiled_f(image, roi_wide, i, j, 0);
            image_cell{i,j} = cutImage;
            check = sum(cutImage(:));
%             disp(['Check = ',num2str(check)])

            % skipping check
            if (check <= size(cutImage,1)*size(cutImage,2)*5) && g == 2
                disp(['Image ',num2str(i),'-',num2str(j),' is skipped because sum of all pixel values is ', num2str(check),'.'])
                break
            end
            
            imwrite(cutImage, [SaveDirPath, imageFilenames{1,g}]);
        end
        
        % skipping check
        if check <= size(cutImage,1)*size(cutImage,2)*5 && g == 2
            continue
        end

        disp(['Cutting process for MPS image ',num2str(i),'-',num2str(j),' done.'])
        
        imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\FerScene\mpsImagesTest\');

        [rows, cols, cameraPoints, projectorPoints, ~] = decodeMPS(imagesCameraDirectoryPath, params);
        
        tformProjector2Camera = fitgeotrans(projectorPoints, cameraPoints, 'projective');
        homography_cell{i,j} = tformProjector2Camera; 
        
        % Resizing the image so it would fit with the estimated homography
        reconstructionResize = imresize(reconstruction,[params.screenHeight, params.screenWidth]);
        RecImageCameraView = imwarp(reconstructionResize, tformProjector2Camera);

        reconstruction_cell{i,j} = RecImageCameraView;
        
        if debug
            figure,
            colormap(gray);
            imagesc(RecImageCameraView)
            pause;
        end
        
        close all
        disp(['Process for MPS image ',num2str(i),'-',num2str(j),' done.'])
        
        toc
    end 
    
    % skipping check
    if check <= size(cutImage,1)*size(cutImage,2)*5 && g == 2
        continue
    end
    
end
    save('result_test.mat','reconstruction_cell','RecImageCameraView','homography_cell','image_cell');
    
% Loading the result if homography is already calculated
else
    load('result_test.mat','reconstruction_cell','RecImageCameraView','homography_cell','image_cell');
end

%% Saving results
% Filling black spots with zeros 
for i = 1:15
    for j = 1:15
        if isempty(reconstruction_cell{i,j})
            reconstruction_cell{i,j} = zeros(size(RecImageCameraView,1),size(RecImageCameraView,2));
%             reconstruction_cell{i,j} = zeros(427,456);
        end
    end
end

% Resizing all images in cell to same size
for i = 1:15
    for j = 1:15
        reconstruction_cell{i,j} = imresize(reconstruction_cell{i,j}, size(RecImageCameraView));
%         reconstruction_cell{i,j} = imresize(reconstruction_cell{i,j}, [427,456]);
    end
end

reconstruction_result = cell2mat(reconstruction_cell);
save('result_matrix.mat', 'reconstruction_result');

imwrite(reconstruction_result,'D:\faks\projekt\script\GitRepository\data\Scenes\result.tiff');

%% Comparing original and warped images
for i = 1:15
    for j = 1:15
        figure(1)
        colormap gray
        subplot(121), imshow(reconstruction_cell{i,j}*1.5);
        axis image
        subplot(122), imshow(cell2{i,j}*1.5);
        axis image
        waitforbuttonpress
    end 
end
