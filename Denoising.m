clear all
close all

% Denoising scipt that calculates median of whole image and replaces all
% pixels with values bigger than treshold * median with median value.
% Histogram calculated for clarification purposes. 

ImagesDirPath = '.\data\imagesCamera\Tiff\';

% variable debug for displaying image histogram and images with and without
% noise, 1 = display, 0 = do not display
debug = 0;

%% Section used for processing whole folder
% dirList = dir(ImagesDirPath);
% isFile = ~[dirList.isdir];
% imageFilenames = {dirList(isFile).name};

%% Section for processing a single picture

    testImg = imread([ImagesDirPath, 'exps_1_100_iso_1600_image_1.tiff']);
    
    if debug == 1
        figure()
        imshow(testImg);
        title('Image with noise');
    
        [vals, edges, bin] = histcounts(testImg);
        figure()
        stem(edges(1:1:end-1), vals);
        title('Histogram');
    end

    medianVal = median(testImg(:));
    thresh = 5;
    hotPixels = double(testImg).*(testImg > (thresh*medianVal));  
    [row, col] = find(hotPixels); 
    testImg = testImg;
    testImg(find(hotPixels)) = medianVal;

    if debug == 1
        figure()
        imshow(testImg);
        title('Image without noise');
    end