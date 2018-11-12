function [DenoisedImage] = Denoising_f(testImg, debug)
 
% Denoising function that calculates median of whole image and replaces all
% pixels with values bigger than treshold * median with median value.
% Histogram calculated for clarification purposes. 

% Variable debug for displaying image histogram and images with and without
% noise, 1 = display, 0 = do not display

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
    thresh = 5;
    hotPixels = double(testImg).*(testImg > (thresh*medianVal));  
    %[row, col] = find(hotPixels);
    testImg(find(hotPixels)) = medianVal;
    DenoisedImage = testImg;
    
    if debug == 1
        figure()
        imshow(DenoisedImage);
        title('Image without noise');
    end
    
end
