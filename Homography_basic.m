SaveDirPath = '.\data\Scenes\FerScene\whiteImagesTest\';
OriginalImage = imread([SaveDirPath, 'white_002.tiff']);

tformP2C = fitgeotrans(projectorPoints, cameraPoints, 'projective');
tformC2P = fitgeotrans(cameraPoints, projectorPoints, 'projective');

[WarpedImage, WarpedView] = imwarp(OriginalImage, tformC2P);
WarpedImage = rescale(WarpedImage, 1080,1920);
% outputView = imref2d(size(OriginalImage));

BackImage = imwarp(WarpedImage, tformP2C);

figure()
colormap gray
subplot(131), imagesc(OriginalImage),title('Original')
axis image
subplot(132), imagesc(WarpedImage),title('Warped')
axis image
subplot(133), imagesc(BackImage),title('Back')
axis image

