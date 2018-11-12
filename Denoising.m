clear all
close all

% Denoising scipt that calculates median of whole image and replaces all
% pixels with values bigger than treshold * median with median value.
% Histogram calculated for clarification purposes. 
%%

fileID = '.\data\imagesCamera\Tiff\exps_1_125_iso_1600_image_1.tiff';
testImg = imread(fileID);

%%
figure()
imshow(testImg);
title('Image with noise');

[vals, edges, bin] = histcounts(testImg);
figure()
stem(edges(1:1:end-1), vals);
title('Histogram');

% index = find(vals == max(vals));
% medianVal = edges(index);

medianVal = median(testImg(:));
thresh = 5;
hotPixels = double(testImg).*(testImg > (thresh*medianVal));  
[row, col] = find(hotPixels); 
testImg = testImg;
testImg(find(hotPixels)) = medianVal;

figure()
imshow(testImg);
title('Image without noise');