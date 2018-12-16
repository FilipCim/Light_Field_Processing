function [outputArg1,outputArg2] = generateDifferentGrayscaleLevelImages(directory, resolution, bitDepth, numberOfLevels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(exist(directory,'dir')~=7)
    mkdir(directory)
end

for i = 1:ceil((2^bitDepth)/numberOfLevels):(2^bitDepth)-1
    grayscaleImage = (i/2^bitDepth).*ones(resolution);

    
    name =  ['grayscaleImage_', num2str(i, '%03d')];
    display(['Writing:     ', name, '.png']);

    imwrite(grayscaleImage, [directory, name, '.png']);
    
end

