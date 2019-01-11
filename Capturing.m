clear all

addpath('.\Lytro Power Tools Beta 1.0.1\lytro-power-tools-1.0.1b0\lpt\bin')

% delay time in seconds
delayTimeSlideShow = 0.5;
delayTimeCamera = 0.2;

% screen ID where to project images
screenID = 2;
screenWidth = 1920;
screenHeight = 1080;

%% Select scene - 1 for White scene, 2 for FER scene

scene = 1;

%%

if scene == 1
    % White scene project images folder paths
    imagesProjectorDirectoryPathBlack = '.\data\Scenes\WhiteScene_all\blackImages\';
    imagesProjectorDirectoryPathWhite = '.\data\Scenes\WhiteScene_all\whiteImages\';
    imagesProjectorDirectoryPathGray = '.\data\Scenes\WhiteScene_all\grayscaleImages\';
    imagesProjectorDirectoryPathMPS = '.\data\Scenes\WhiteScene_all\MPSImages\';
    
    % White scene RAW images folder paths (save folders)
    imagesCameraDirectoryPathRawBlack = '.\data\Scenes\WhiteScene_all\blackImagesRaw\';
    imagesCameraDirectoryPathRawWhite = '.\data\Scenes\WhiteScene_all\whiteImagesRaw\';
    imagesCameraDirectoryPathRawGray = '.\data\Scenes\WhiteScene_all\grayscaleImagesRaw\';
    imagesCameraDirectoryPathRawMPS = '.\data\Scenes\WhiteScene_all\MPSImagesRaw\';
    
    dirList = dir(imagesProjectorDirectoryPathGray);
    isFile = ~[dirList.isdir];
    imageFilenames3 = {dirList(isFile).name};
end

if scene == 2
    % FER scene project images folder paths
    imagesProjectorDirectoryPathBlack = '.\data\Scenes\WhiteScene_all\blackImages\';
    imagesProjectorDirectoryPathWhite = '.\data\Scenes\WhiteScene_all\whiteImages\';
    imagesProjectorDirectoryPathHadamard = '.\data\Scenes\FerScene\hadamardImages\';
    imagesProjectorDirectoryPathMPS = '.\data\Scenes\WhiteScene_all\MPSImages\';
    
    % FER scene RAW images folder paths (save folders)
    imagesCameraDirectoryPathRawBlack = '.\data\Scenes\FerScene\blackImagesRaw\';
    imagesCameraDirectoryPathRawWhite = '.\data\Scenes\FerScene\whiteImagesRaw\';
    imagesCameraDirectoryPathRawHadamard = '.\data\Scenes\FerScene\hadamardImagesRaw\';
    imagesCameraDirectoryPathRawMPS = '.\data\Scenes\FerScene\MPSImagesRaw\';
    
    dirList = dir(imagesProjectorDirectoryPathHadamard);
    isFile = ~[dirList.isdir];
    imageFilenames4 = {dirList(isFile).name};
end

dirList = dir(imagesProjectorDirectoryPathBlack);
isFile = ~[dirList.isdir];
imageFilenames1 = {dirList(isFile).name};

dirList = dir(imagesProjectorDirectoryPathWhite);
isFile = ~[dirList.isdir];
imageFilenames2 = {dirList(isFile).name};

dirList = dir(imagesProjectorDirectoryPathMPS);
isFile = ~[dirList.isdir];
imageFilenames5 = {dirList(isFile).name};

%% Setting up Lytro parameters

% Locking exposition
exposition =60;
system('python cameracontrols.py exposureMode -s manual')
system(['python cameracontrols.py shutterSpeed -s 1/',num2str(exposition)]);

% Locking focus
% First measurement was done by manually setting focus to 'auto' to let
% camera focus in the scene and then locking it programmatically. 

% Meanwhile, we saw that there is a option to programmatically set focus to
% 'auto', but I'm not sure if that will keep changing focus after locking
% it. That should be tested.

% system('python cameracontrols.py focusMode -s auto');

system('python cameracontrols.py focusMode -s manual');
system('python cameracontrols.py focusLock -s 1');

% Locking ISO value
isovar=80;
system(['python cameracontrols.py iso -s ', num2str(isovar)]);

system('python cameracontrols.py focusLock -s --1');

%% blackImages
tic
for i=1:8
    image = imread([imagesProjectorDirectoryPathBlack, imageFilenames1{i}]);
    fullscreen(image, screenID);
    
    if(i==1)
        pause(1.3*delayTimeSlideShow)
        else
        pause(delayTimeSlideShow)
    end
    display(['Slika crna broj ',  num2str(i)]);
   system('python mycontrol.py')
   
end
 closescreen;
 display('Skidanje...');
 system(['python mydownload.py download-images --nImages ',num2str(8),' --path ', imagesCameraDirectoryPathRawBlack])
toc
 %% whiteImages
 tic
 for i=1:8
    image = imread([imagesProjectorDirectoryPathWhite, imageFilenames2{i}]);
    fullscreen(image, screenID);
    
    if(i==1)
        pause(1.3*delayTimeSlideShow)
        else
        pause(delayTimeSlideShow)
    end
    display(['Slika bijela broj ',  num2str(i)]);
   system('python mycontrol.py')
   
end
 closescreen;
  display('Skidanje...');
 system(['python mydownload.py download-images --nImages ',num2str(8),' --path ', imagesCameraDirectoryPathRawWhite])
 toc
 
 %% grayscaleImages
 tic
 
 for i=1:8:256
    image = imread([imagesProjectorDirectoryPathGray, imageFilenames3{i}]);
    fullscreen(image, screenID);
    
    if(i==1)
        pause(1.3*delayTimeSlideShow)
        else
        pause(delayTimeSlideShow)
    end
    display(['Slika siva broj ',  num2str(i)]);
   system('python mycontrol.py')
   
end
 closescreen;
 
 disp('Skidanje...')
 system(['python mydownload.py download-images --nImages ',num2str(256/8),' --path ', imagesCameraDirectoryPathRawGray])
 toc
 %% hadamardMeasurement
 tic
 for i=1:256
    image = imread([imagesProjectorDirectoryPathHadamard, imageFilenames4{i}]);
    fullscreen(image, screenID);
    
    if(i==1)
        pause(1.3*delayTimeSlideShow)
        else
        pause(delayTimeSlideShow)
    end
    display(['Slika hadamard broj ',  num2str(i)]);
   system('python mycontrol.py')
   
end
 closescreen;
  display('Skidanje');
 system(['python mydownload.py download-images --nImages ',num2str(256),' --path ', imagesCameraDirectoryPathRawHadamard])
 toc
 
 
 %% mps
 tic
  for i=1:50
    image = imread([imagesProjectorDirectoryPathMPS, imageFilenames5{i}]);
    fullscreen(image, screenID);
    
    if(i==1)
        pause(1.3*delayTimeSlideShow)
        else
        pause(delayTimeSlideShow)
    end
    display(['Slika mps broj ',  num2str(i)]);
   system('python mycontrol.py')
   
  end
 closescreen;
  display('Skidanje');
 system(['python mydownload.py download-images --nImages ',num2str(50),' --path ', imagesCameraDirectoryPathRawMPS])
 toc