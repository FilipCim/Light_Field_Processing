clear all
close all

% cuts ROI out of the images and concatenates them into a single tiled TIF
% image

debug = 1;

%% Section for concatenating images together into a TIFF image

% Don't have necessarry metadata from camera so the testing is done on one
% of the images from the sample data

% "fake data"
ImagesDirPathLF = '.\data\Scenes\whiteScene_all\whiteImagesRaw\';
LFP = 'IMG_0904.LFR';

% White image decoding
WhiteImagesDirPath = '.\data\WhiteImages\Cameras\';
LFUtilProcessWhiteImages(WhiteImagesDirPath);
    
% LFR decoding
DecodeOptions.WhiteImageDatabasePath = [WhiteImagesDirPath, 'WhiteImageDatabase.mat'];
LF = LFLytroDecodeImage([ImagesDirPathLF, LFP], DecodeOptions);

%% using one image to determine ROI
img = squeeze(LF(7,7,:,:,1));

roi = ROI_f(img, debug);
roi_wide = round(roi.*1.2);

%% cutting all images and concatenating them into one image

for i=1:15
    for j=1:15
        % differences in real data - squeeze(LF.Raw())./2^10
        im_cell{i,j} = (squeeze(LF(j,i,:,:,1)));
        % roi used instead of roi_wide because of the position of ROI in
        % image, use roi_wide on real data
        % ROI is widened by 20% because of the pixel offset in Lytro camera
        im_cell_conc{i,j} = im_cell{i,j}(floor(roi(2)) : floor(roi(2))...
             + roi(4), floor(roi(1)) : floor(roi(1)) + roi(3));
    end
end

%%
figure()
title('Original tiled image');
lf_tiled = cell2mat(im_cell);
imshow(lf_tiled);

figure()
title('Croped tiled image');
lf_tiled_conc = cell2mat(im_cell_conc);
imshow(lf_tiled_conc);


SaveDirPath = '.\data\imagesCamera\Tiled_tiff\';
imwrite(double(lf_tiled_conc), [SaveDirPath, 'image_tiled.tiff']);
