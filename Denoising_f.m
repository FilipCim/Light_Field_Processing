function [DenoisedImage, row, col] = Denoising_f(testImg, debug, row, col)
 
% Denoising function that calculates median of whole image and replaces all
% pixels with values bigger than treshold * median with median value.
% Histogram calculated for clarification purposes.

% ------------- add [row, col] inputs so that hot pixel finding is not done
% every time------------

% variable black determins whether the hot pixels are going to be calculated
% or not - 1 = calculate, 0 = don't calculate.

% Variable debug for displaying image histogram and images with and without
% noise, 1 = display, 0 = do not display

if nargin == 1
    debug = 0;
end

%% Section for processing a single picture
    
    if debug == 1
        figure()
        imshow(testImg);
        title('Image with noise');
    
        [vals, edges] = histcounts(testImg);
        figure()
        stem(edges(1:1:end-1), vals);
        title('Histogram');
    end

    medianVal = median(testImg(:));
    
    if nargin <= 2
        thresh = 5;
        hotPixels = double(testImg).*(testImg > (thresh*medianVal)); 
        [row, col] = find(hotPixels);
    end
    
    testImg(row,col) = medianVal;
    DenoisedImage = testImg;
    
    if debug == 1
        figure()
        imshow(DenoisedImage);
        title('Image without noise');
    end
    
end
