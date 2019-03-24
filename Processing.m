clear all
close all
clc

utilities_functions = './utilities';
mps_functions = './utilities/mps';
addpath(utilities_functions, mps_functions); 

% SeDuMi path, used for initial testing (already precompiled)
addpath('C:\development\SeDuMi_1_3\');
addpath('C:\development\spams-matlab-v2.6');

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

% variable for drawing some of the additional plots. 1 = draw, 0 = don't
debug = 1;
%% Calculating ROI 
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
        figure(1122)
        imshow(img);
        hold on;
        rectangle('Position',...
            [roi_wide(1),roi_wide(2),roi_wide(3),roi_wide(4)],...
            'EdgeColor','r','LineWidth',2 )
        drawnow
        hold off;
    end

disp('Found ROI...')

%%
% Prealocating cell for hadamard images 
hadamard_cell = cell(15,15,256);

% Prealocating cell for reconstructed images
reconstruction_cell = cell(15,15);

%% Selecting ROI on the image
imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\mpsImagesTest\');

whiteImageScene = imread([imagesCameraDirectoryPath, 'white', '.tiff']);
whiteImageScene = double(whiteImageScene)./double(max(whiteImageScene(:)));

firstPhaseImage = ((imread(fullfile(currentDirectory,... 
    '\data\Scenes\FerScene\hadamardImagesTiled\measurement_mask_1_01.tiff'))));
firstPhaseImage = double(firstPhaseImage);

firstPhaseImage = CutTiled_f(firstPhaseImage, roi_wide, 7, 7, 0);
 
    figure(111),
    imshowpair(whiteImageScene, firstPhaseImage)
    figure(111)

    [xi, yi] = getpts(gcf);
    close gcf;

disp('ROI selected.')

%%
% Two loops for iterating through rows and columns of tiled images. i =
% rows, j = columns

% skipping black images
imageCounti = []

%for i = 1:15
% for testing
for i = 1:1
    %for j = 1:15
    % for testing
    for j = 9:9

% Whole processing is done on every tile of tiled image (evey pixel view)
tic

%% CALCULATE PER PIXEL CORRESPONDENCES USING MPS FOR CALIBRATION
% Cutting tiled MPS image and saving it into a folder
TiledImagePath = '.\data\Scenes\WhiteScene_all\mpsImagesTiled\';
SaveDirPath = '.\data\Scenes\WhiteScene_all\mpsImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all MPS images and for one pixel view (one tile in tiled image)
for g = 1:size(imageFilenames, 2)
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    cutImage = CutTiled_f(image, roi_wide, i, j, 0);
    imwrite(cutImage, [SaveDirPath, imageFilenames{1,g}]);
end

disp('Cutting process for MPS images done (White scene).')

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
% Cutting white images for devignetting
TiledImagePath = '.\data\Scenes\WhiteScene_all\whiteImagesTiled\';
SaveDirPath = '.\data\Scenes\WhiteScene_all\whiteImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all white images and for one pixel view (one tile in tiled image)
for g = 1:size(imageFilenames, 2)
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    cutImage = CutTiled_f(image, roi_wide, i, j, 0);
% imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\whiteImagesTest\');
% 
% dirList = dir(imagesCameraDirectoryPath);
% isFile = ~[dirList.isdir];
% imageFilenames = sort({dirList(isFile).name});

    imagesCamera(:,:,g) = cutImage;
%     figure(100)
%     imagesc(imagesCamera(:,:,i));
%     waitforbuttonpress
%     axis image
end

whiteImageCamera= mean(imagesCamera,3);
whiteImageProjectorWarped = imwarp(ones([params.screenHeight, params.screenWidth]), tformProjector2Camera, 'OutputView',outputView);

whiteImageCamera = double(whiteImageCamera)./double(max(whiteImageCamera(:)));

vignettingFunctionEstimation = double(whiteImageProjectorWarped)./whiteImageCamera;

figure(),
imagesc(vignettingFunctionEstimation)
title('Estimated vignetting function')

%%
% Importing white image...
imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\mpsImagesTest\');

whiteImageScene = imread([imagesCameraDirectoryPath, 'white', '.tiff']);
whiteImageScene = double(whiteImageScene)./double(max(whiteImageScene(:)));
% % 
% % figure()
% % imagesc(whiteImageScene.*vignettingFunctionEstimation)

%% ESTIMATE NOISE
% Cutting black images for devignetting
TiledImagePath = '.\data\Scenes\WhiteScene_all\blackImagesTiled\';
SaveDirPath = '.\data\Scenes\WhiteScene_all\blackImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all black images and for one pixel view (one tile in tiled image)
for g = 1:size(imageFilenames, 2)
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    cutImage = CutTiled_f(image, roi_wide, i, j, 0);

% imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\whiteScene_all\blackImagesTest\');
% 
% dirList = dir(imagesCameraDirectoryPath);
% isFile = ~[dirList.isdir];
% imageFilenames = sort({dirList(isFile).name});

    imagesCamera(:,:,g) = cutImage;
    imagesCameraWarped(:,:,g) = imagesCamera(:,:,g);
    
end

% noiseEstimation = mean(imagesCameraWarped, 3).*vignettingFunctionEstimation;
noiseEstimation = mean(imagesCameraWarped, 3);

figure(),
imagesc(noiseEstimation)
title('Estimated noise')


%% CALCULATE PER PIXEL CORRESPONDENCES USING MPS FOR SCENE

TiledImagePath = '.\data\Scenes\FerScene\mpsImagesTiled\';
SaveDirPath = '.\data\Scenes\FerScene\mpsImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all MPS images and for one pixel view (one tile in tiled image)
for g = 1:size(imageFilenames, 2)
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    cutImage = CutTiled_f(image, roi_wide, i, j, 0);
    imwrite(cutImage, [SaveDirPath, imageFilenames{1,g}]);
end

disp('Cutting process for MPS images done (FER scene).')


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
% firstPhaseImage = ((imread(fullfile(currentDirectory,... 
%     '\data\Scenes\FerScene\hadamardImagesTiled\measurement_mask_1_01.tiff'))));
% firstPhaseImage = double(firstPhaseImage);
% 
% firstPhaseImage = CutTiled_f(firstPhaseImage, roi_wide, i, j, 0);
% 
% 
%     figure(111),
%     imshowpair(whiteImageScene, firstPhaseImage)
%     figure(111)
% 
%     [xi, yi] = getpts(gcf);
%     close gcf;


% ROI selection done but manually inputing two points on the image. First
% one is upper left. Second one is lower right. Last point needs to be
% double-click, shit-click or right click to be able to add the last point.

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
    for h = 1:size(idx,2)
        pixValue = repmat(impixel(phases{phaseNo}, centroidCoords(h,1),centroidCoords(h,2)), size(phasesWatershed{phaseNo}(sWater(idx(h)).PixelIdxList)));
        
        phases{phaseNo}(sWater(idx(h)).PixelIdxList) = pixValue(:,1);
    end
    
    
end




%%

if debug
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
end

%% Measuring
% Cutting hadamard images and putting them in a cell
TiledImagePath = '.\data\Scenes\FerScene\hadamardImagesTiled\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

for g = 1:size(imageFilenames, 2)           
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    hadamard_cell{i,j,g} = CutTiled_f(image, roi_wide, i, j, 0);
end

disp('Cutting hadamard images done')

%%
% imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\FerScene\hadamardImagesTest\');

% Processing phase-by-phase
% Here we iterate every phase for every measurement (in the hadamard 
% folder)

counter = 0;

for phaseNumber = 1:4
    for measurementNumber = 1:64
        display(['Processing:     (',num2str((phaseNumber-1)*64+measurementNumber),'/',num2str(4*64),')', imagesCameraDirectoryPath, 'measurement_mask_', num2str(phaseNumber), '_', num2str(measurementNumber, '%02d'), '.tiff'])
        % counter for reading images from a cell
        counter = counter+1;
        image = hadamard_cell{i,j,counter};
        image = double(image);
        
        % Noise reduction
%         image = double(image).*vignettingFunctionEstimation;
        image = image - noiseEstimation;    
        
        %tic
        
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
                        measurements(row-startRow+1, col-startCol+1, measurementNumber)= mean2((image(mask)));
% 
%                         figure(101)
%                         imshowpair(image, mask)
%                         drawnow
                        
                    end
                end
                
   
        end
                
        %toc
    end
end

display('Measurement preprocessing ended!')

%% Check measurments
if debug
figure()
colormap gray
imagesc(measurements(:,:,22));
title('One of the measurements');
end

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

ind = 2:64;

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
reconstructionOptions.Package = 'sedumi';
% for cvx reconstruction package user has choice of algorithms(sedumi, sdpt3)
% reconstructionOptions.Algorithm = 'sedumi';
reconstructionOptions.lambda = 0;

reconstruction=[];



for i = 1:size(measurements, 1)
    for j = 1:size(measurements, 2)
        
        subimage_estimation = l1Optimization(y(i, j, :), psi_inv, theta, reconstructionOptions);
        
        subimage_estimation = reshape(subimage_estimation, block_size, block_size);
        reconstruction((i-1)*block_size+1:(i-1)*block_size+block_size,(j-1)*block_size+1:(j-1)*block_size+block_size) = subimage_estimation;
        
    end
    
    if debug
    figure(101)
    imagesc(reconstruction), axis image, colormap gray
    drawnow
    end
end

reconstruction_cell{i,j} = reconstruction;

if debug
    figure(1011)
    imagesc(reconstruction_cell{i,j}), axis image, colormap gray
    drawnow    
end
%% CALCULATE CALIBRATION PATCHES


coeffOfVar = std(measurements./repmat(mean(measurements,3), [1,1,64]),0,3);
coeffOfVar = coeffOfVar./max(coeffOfVar(:));

% figure(), imagesc(coeffOfVar)
% imagesc(coeffOfVar<0.975*median(coeffOfVar(:)))

display(['Reconstructed image ',num2str(i),' // ', num2str(j)])
toc

% Ending two main loops that iterate through rows and columns of tiled
% image.
end
end


%% Saving image
% SaveDirPath = '.\data\Scenes\Results';
% reconstructedImage = cell2mat(reconstruction_cell);
% 
% imwrite(double(lf_tiled_conc), [SaveDirPath, 'Recontructed_Image.tiff']);
% 
% disp('Done!')


