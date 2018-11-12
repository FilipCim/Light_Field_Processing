clearvars
close all
clc

% add utilities path to MATLAB path
% utilities contain commonly used scipts
run('C:\Users\Nicol\Downloads\FER\8. SEM\diplomski seminar\LF Matlab\LF Toolbox\LFMatlabPathSetup.m');
addpath('./utilities')
addpath('./LF_Pictures')
addpath('./data/imagesCamera/Raw/')
addpath('C:\PROJEKT\lytro\Lytro Power Tools Beta 1.0.1\lytro-power-tools-1.0.1b0\lpt\bin')


%system('python cameracontrols.py exposureMode -g --get');
%%

% get all image filenames from a directory
% these are to be projected
%imagesProjectorDirectoryPath = './data/Playing_Cards/';
imagesProjectorDirectoryPath = './data/mpsImages/';
imagesCameraDirectoryPathtest='./data/imagesCamera';
imagesCameraDirectoryPathTiff = './data/imagesCamera/Tiff/';
imagesCameraDirectoryPathRaw= './data/imagesCamera/Raw';

if(exist(imagesCameraDirectoryPathtest, 'dir')~=7)
    mkdir imagesCameraDirectoryPathtest
end   
if(exist(imagesCameraDirectoryPathTiff, 'dir')~=7)
    mkdir './data/imagesCamera/Tiff'
end
if(exist(imagesCameraDirectoryPathRaw, 'dir')~=7)
    mkdir './data/imagesCamera/Raw'
end
    
dirList = dir(imagesProjectorDirectoryPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


% delay time in seconds
delayTimeSlideShow = 0.5;
delayTimeCamera = 0.2;

% screen ID where to project images
screenID = 1;
screenWidth = 1920;
screenHeight = 1080;


% % these properties are related to GigE vison camera
% % you'll need to change these to correspond to Lytro properties
% cam = gigecam;
% cam.ExposureTime = 1666.6600341796875;
% cam.AcquisitionFrameRate = 5;
% cam.AcquisitionFrameRateAbs = 5;
% cam.GainAutoBalance = 'off';
% cam.BlackLevelAutoBalance = 'off'

system('python cameracontrols.py exposureMode -s --manual 1/64 ')


%%

% sets black screen on the defined screenID
%fullscreen(zeros([screenWidth, screenHeight]), screenID)

% this loop shows all the images  in the folder on the defined screenID
% and captures them on the camera

%for i = 1:size(imageFilenames, 2)
for i=1:50 

    tic
    
    % load image
    image = imread([imagesProjectorDirectoryPath, imageFilenames{i}]);
    image = imresize(image, [1080, 1920]);
    
    % set it fullscreen
    fullscreen(image, screenID);
    
    % add a delay so the image stabilizes on the display - customizable
    if(i==1)
        pause(1.3*delayTimeSlideShow)
    else
        pause(delayTimeSlideShow)
    end
    
    tic
   system('python mycontrol.py')
    toc
    
end

closescreen

tic
system('python mydownload.py download-images --nImages 50')
toc

  %654 camera bin

%% raspakiravanje
tic
dirList = dir(imagesCameraDirectoryPathRaw);
isFile = ~[dirList.isdir];
imageFilenames2 = {dirList(isFile).name};

for i=1:50
raw_image=LFReadLFP(imageFilenames2{i});
proba=raw_image.RawImg;
proba = double(proba)./ 2^10;

filename=imageFilenames{i}
[path, name, ext]=fileparts(filename)
imwrite(proba, [imagesCameraDirectoryPathTiff, name,'.tiff']);
end
toc










