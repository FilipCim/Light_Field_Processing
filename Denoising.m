clear all
close all
clc

%%

testImg = load(

[vals, edges, bin] = histcounts(testImg);
thresh = 5; 
hotPixels = double(testImg).*(testImg>(thresh*medianVal));  
[row, col] = find(hotPixels); 
testImg(find(hotPixels))=medianVal;
