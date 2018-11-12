clearvars
close all
clc

t = tic;
%%
% add utilities path to MATLAB path
% utilities contain commonly used scipts
run('C:\Users\Nicol\Downloads\FER\8. SEM\diplomski seminar\LF Matlab\LF Toolbox\LFMatlabPathSetup.m');
addpath('./utilities')
addpath('./LF_Pictures')
addpath('./data/imagesCamera/Raw/')

% get all image filenames from a directory
% these are to be projected
%imagesProjectorDirectoryPath = './data/Playing_Cards/';
%imagesProjectorDirectoryPath = './data/mpsImages/';
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

%%
isovar=1600;
system(['python cameracontrols.py iso -s ', num2str(isovar)]);
system('python cameracontrols.py exposureMode -s manual')
exps=[125 100 80 60 50 40 30 25 20 15 12 10 8 6.4 5 4];

for i=1:length(exps)
         system(['python cameracontrols.py shutterSpeed -s 1/',num2str(exps(i))]);
        for i=1:4
        system('python mycontrol.py')
        end
end 

system(['python mydownload.py download-images --nImages ',num2str(length(exps)*4),' --path ', imagesCameraDirectoryPathRaw])

Time1 = toc(t)


%% raspakiravanje
dirList = dir(imagesCameraDirectoryPathRaw);
isFile = ~[dirList.isdir];
imageFilenames2 = {dirList(isFile).name};

tic
for j=1:length(exps)
    
    for i=1:4    
    raw_image=LFReadLFP(imageFilenames2{i});
    proba=raw_image.RawImg;
    proba = double(proba)./ 2^10;
    imwrite(proba, [imagesCameraDirectoryPathTiff, 'exps_1_',num2str(exps(j)),'_iso_',num2str(isovar),'_image_',num2str(i),'.tiff']);
    end
    
end
toc
%Time2 = toc(t)