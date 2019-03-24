close all

% Denoising scipt that calculates median of whole image and replaces all
% pixels with values bigger than treshold * median with median value.
% Histogram calculated for clarification purposes. 

%%
% ImagesDirPath = '.\data\Scenes\whiteScene_all\blackImagesRaw\';
% LFR_image=LFReadLFP([ImagesDirPath, 'IMG_0896.LFR']);
% raw_image=LFR_image.RawImg;
% raw_image = double(raw_image)./ 2^10;
% 
% imagesCameraDirectoryPathTiff = '.\data\imagesCamera\Tiff\';
% imwrite(raw_image, [imagesCameraDirectoryPathTiff, 'black_test','.tiff']);
%%
ImagesDirPath = '.\data\imagesCamera\Tiff\';

% variable debug for displaying image histogram and images with and without
% noise, 1 = display, 0 = do not display
debug = 1;
black = 1;

%% Section for processing a single picture

    testImg = imread([ImagesDirPath, 'black.tiff']);
    
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
    
    if black == 1
        thresh = 5;
        hotPixels = double(testImg).*(testImg > (thresh*medianVal)); 
        [row, col] = find(hotPixels);
    end
    
    testImg(row,col) = medianVal;
    DenoisedImage = testImg;

    if debug == 1
        figure()
        imshow(testImg);
        title('Image without noise');
    end