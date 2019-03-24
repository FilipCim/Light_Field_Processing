close all

% This script unpacks all of the LFR files and turnes them into tiled TIFF
% images saving them in into appropriate folders. MPS images are turned
% into tiled images and every tiled look (cell part) is saved into their
% folder (the whole series of MPS images are saved in one folder).

% The same structure will be used for unpacking every folder of images,
% only the source, destination and names change. This is not the case only
% for the special look-by-look unpacking of the MPS images

% Unpacking will be done first for the white scene, and then for the FER
% scene

%% White scene
%% Finding ROI and White Images Decoding
% White image decoding
WhiteImagesDirPath = '.\data\WhiteImages\Cameras\';
LFUtilProcessWhiteImages(WhiteImagesDirPath);

DecodeOptions.WhiteImageDatabasePath = [WhiteImagesDirPath, 'WhiteImageDatabase.mat'];

%% Finding ROI
% WhiteDirPath = '.\data\Scenes\whiteScene_all\whiteImagesRaw\';
% WhiteImage = 'IMG_1759.LFR';
%  
% LF = LFLytroDecodeImage([WhiteDirPath, WhiteImage], DecodeOptions);
% 
% img = squeeze(LF(7,7,:,:,1));
% save('.\data\Scenes\ROI_image2.mat', 'img')
%%
% Loading preunpacked Light Field of white image for ROI finding
ROI_image = load('.\data\Scenes\ROI_image2Kosa.mat');
img = ROI_image.img;
roi = ROI_f(img, 1);

% Increasing ROI for 20% and checking the size of the ROI so it would not
% extend the size of the image

% 10% of the width (x dim)
roi_scale(1) = roi(3)*0.1;
% 10% of the length (y dim)
roi_scale(2) = roi(4)*0.1;

roi_wide(1) = round(roi(1) - roi_scale(1));
if roi_wide(1) < 0
    roi_wide(1) = 1;
end

roi_wide(2) = round(roi(2) - roi_scale(1));
if roi_wide(2) < 0
    roi_wide(2) = 1;
end

roi_wide(3) = round(roi(3) + 2*roi_scale(1));
if (roi_wide(1)+roi_wide(3)) > size(img, 2)
    roi_wide(3) = size(img, 1) - roi_wide(1);
end

roi_wide(4) = round(roi(4) + 2*roi_scale(2));
if (roi_wide(4)+roi_wide(2)) > size(img, 1)
    roi_wide(4) = size(img, 1) - roi_wide(2);
end

debug = 1;
if debug
        figure()
        imshow(img);
        hold on;
        rectangle('Position', [roi_wide(1),roi_wide(2),roi_wide(3),roi_wide(4)],'EdgeColor','r','LineWidth',2 )
    end

disp('Found ROI...')

%% 
disp('.')
disp('.')
disp('.')
disp('Decoding White scene')
disp('.')
disp('.')
disp('.')

%% Black images
% some slight modification are done in the naming aspect. Names dedicated
% are black_x with x being the number of the image in folder

SourceDirPath = '.\data\Scenes\whiteScene_all\blackImagesRaw\';
SaveDirPath = '.\data\Scenes\whiteScene_all\blackImagesTiled2\';
% NamesDirPath = '.\data\Scenes\whiteScene_all\blackImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


% dirList2 = dir(NamesDirPath);
% isFile2 = ~[dirList2.isdir];
% imageFilenames2 = {dirList2(isFile2).name};

% decoding black images
disp('Decoding black images')
%for g=1:size(imageFilenames, 2)
for g=1:1    
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

% filename = imageFilenames2{g};
% [path, name, ext]=fileparts(filename);
name = ['black_', sprintf('%03d',g)];
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end

%% Grayscale images
SourceDirPath = '.\data\Scenes\whiteScene_all\grayscaleImagesRaw\';
SaveDirPath = '.\data\Scenes\whiteScene_all\grayscaleImagesTiled\';
NamesDirPath = '.\data\Scenes\whiteScene_all\grayscaleImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


dirList2 = dir(NamesDirPath);
isFile2 = ~[dirList2.isdir];
imageFilenames2 = {dirList2(isFile2).name};
imageFilenames2 = imageFilenames2(1:8:end);

% Decoding grayscale images
disp('Decoding grayscale images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

% multiplying the loop index by 8 beacause the images were taken with the
% offset of 8
% this way every image has an appropriate name

filename = imageFilenames2{1,g};
[path, name, ext]=fileparts(filename);
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end


%% MPS images tiled
% this part of the script produces a tiled image of MPS images
% not sure if this is needed (but I have it here so I'll just leave it)
SourceDirPath = '.\data\Scenes\whiteScene_all\mpsImagesRaw\';
SaveDirPath = '.\data\Scenes\whiteScene_all\mpsImagesTiled2\';
NamesDirPath = '.\data\Scenes\whiteScene_all\mpsImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


dirList2 = dir(NamesDirPath);
isFile2 = ~[dirList2.isdir];
imageFilenames2 = {dirList2(isFile2).name};

% Decoding MPS images
disp('Decoding MPS images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

filename = imageFilenames2{g};
[path, name, ext]=fileparts(filename);
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end

%% White images
% some slight modification are done in the naming aspect. Names dedicated
% are white_x with x being the number of the image in folder

SourceDirPath = '.\data\Scenes\whiteScene_all\whiteImagesRaw\';
SaveDirPath = '.\data\Scenes\whiteScene_all\whiteImagesTiled2\';
% NamesDirPath = '.\data\Scenes\whiteScene_all\whiteImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


% dirList2 = dir(NamesDirPath);
% isFile2 = ~[dirList2.isdir];
% imageFilenames2 = {dirList2(isFile2).name};

% Decoding whtie images
disp('Decoding white images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

% filename = imageFilenames2{g};
% [path, name, ext]=fileparts(filename);
name = ['white_', sprintf('%03d',g)];
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end

%% ================================================================================================================
%% ================================================================================================================
%% FER scene
disp('.')
disp('.')
disp('.')
disp('Decoding FER scene')
disp('.')
disp('.')
disp('.')

%% Black images
% some slight modification are done in the naming aspect. Names dedicated
% are black_x with x being the number of the image in folder

SourceDirPath = '.\data\Scenes\FerScene\blackImagesRaw\';
SaveDirPath = '.\data\Scenes\FerScene\blackImagesTiled\';
% NamesDirPath = '.\data\Scenes\whiteScene_all\blackImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


% dirList2 = dir(NamesDirPath);
% isFile2 = ~[dirList2.isdir];
% imageFilenames2 = {dirList2(isFile2).name};

% decoding black images
disp('Decoding black images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

% filename = imageFilenames2{g};
% [path, name, ext]=fileparts(filename);
name = ['black_', sprintf('%03d',g)];
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end

%% Hadamard images 
% this part of the script decodes hadamard images used for measuring.

SourceDirPath = '.\data\Scenes\FerScene\hadamardImagesRaw\';
SaveDirPath = '.\data\Scenes\FerScene\hadamardImagesTiled\';
NamesDirPath = '.\data\Scenes\FerScene\hadamardImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


dirList2 = dir(NamesDirPath);
isFile2 = ~[dirList2.isdir];
imageFilenames2 = {dirList2(isFile2).name};

% Decoding hadamard images
disp('Decoding hadamard images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

filename = imageFilenames2{g};
[path, name, ext]=fileparts(filename);
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end

%% MPS images tiled
% this part of the script produces a tiled image of MPS images
% not sure if this is needed (but I have it here so I'll just leave it)
SourceDirPath = '.\data\Scenes\FerScene\mpsImagesRaw\';
SaveDirPath = '.\data\Scenes\FerScene\mpsImagesTiled\';
NamesDirPath = '.\data\Scenes\whiteScene_all\mpsImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


dirList2 = dir(NamesDirPath);
isFile2 = ~[dirList2.isdir];
imageFilenames2 = {dirList2(isFile2).name};

% Decoding MPS images
disp('Decoding MPS images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

filename = imageFilenames2{g};
[path, name, ext]=fileparts(filename);
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end


%% White images
% some slight modification are done in the naming aspect. Names dedicated
% are white_x with x being the number of the image in folder

SourceDirPath = '.\data\Scenes\FerScene\whiteImagesRaw\';
SaveDirPath = '.\data\Scenes\FerScene\whiteImagesTiled\';
% NamesDirPath = '.\data\Scenes\whiteScene_all\whiteImages\';

dirList = dir(SourceDirPath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};


% dirList2 = dir(NamesDirPath);
% isFile2 = ~[dirList2.isdir];
% imageFilenames2 = {dirList2(isFile2).name};

% Decoding whtie images
disp('Decoding white images')
for g=1:size(imageFilenames, 2)
LF = LFLytroDecodeImage([SourceDirPath, imageFilenames{g}], DecodeOptions);

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi_wide(2)) : floor(roi_wide(2))...
             + roi_wide(4), floor(roi_wide(1)) : floor(roi_wide(1)) + roi_wide(3));
    end
end

disp('Made a tiled image...')

lf_tiled_conc = cell2mat(im_cell_conc);

% filename = imageFilenames2{g};
% [path, name, ext]=fileparts(filename);
name = ['white_', sprintf('%03d',g)];
imwrite(double(lf_tiled_conc), [SaveDirPath, name,'.tiff']);

disp(['Finished image ', name, '.tiff'])
disp('=========================Sheep=====================================')
end