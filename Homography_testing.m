%% Scaling reconstruction
reconstruction = mat2gray(reconstruction);

%% INVERSE HOMOGRAPHY
% All of measurements and reconstruction are made in the Projector plane.
% With homography (inverse one) plane of viewing is shifted to the Camera
% plane.
% This should be done for every PixelView image in the tiled format.
% Reconstruction is the same, only thing that is different is homography,
% so MPS should be calculated for every image, but reconstruction shouldn't
debug = 0;
homograph = 1;

% Prealocating cells 
reconstruction_cell = cell(15,15);
homography_cell = cell(15,15);
image_cell = cell(15,15);

% Variable homograph dictates whether or not the (inverse) homography is
% calculated
if homograph
    tic
% Calculating ROI 
% Loading unpacked Light Field (here so the unpacking wouldn't be done
% every time ROI is calculated)
ROI_image = load('.\data\Scenes\ROI_image2Kosa.mat');
img = ROI_image.img;
roi = ROI_f(img, debug);

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

if debug
        figure(112)
        imshow(img);
        hold on;
        rectangle('Position',...
            [roi_wide(1),roi_wide(2),roi_wide(3),roi_wide(4)],...
            'EdgeColor','r','LineWidth',2 )
        drawnow
        hold off;
end

disp('Found ROI...')

% CALCULATE PER PIXEL CORRESPONDENCES USING MPS FOR SCENE
for i = 1:15
    for j = 1:15
        tic
        
        TiledImagePath = '.\data\Scenes\FerScene\mpsImagesTiled\';
        SaveDirPath = '.\data\Scenes\FerScene\mpsImagesTest\';

        dirList = dir(TiledImagePath);
        isFile = ~[dirList.isdir];
        imageFilenames = {dirList(isFile).name};

        % For all MPS images and for one pixel view (one tile in tiled image)
        for g = 1:size(imageFilenames, 2)
            image = imread([TiledImagePath, imageFilenames{1,g}]);
            cutImage = CutTiled_f(image, roi_wide, i, j, 0);
            image_cell{i,j} = cutImage;
            check = sum(cutImage(:));
%             disp(['Check = ',num2str(check)])

            % skipping check
            if (check <= size(cutImage,1)*size(cutImage,2)*5) && g == 2
                disp(['Image ',num2str(i),'-',num2str(j),' is skipped because sum of all pixel values is ', num2str(check),'.'])
                break
            end
            
            imwrite(cutImage, [SaveDirPath, imageFilenames{1,g}]);
        end
        
        % skipping check
        if check <= size(cutImage,1)*size(cutImage,2)*5 && g == 2
            continue
        end

        disp(['Cutting process for MPS image ',num2str(i),'-',num2str(j),' done.'])
        
        imagesCameraDirectoryPath = fullfile(currentDirectory,'\data\Scenes\FerScene\mpsImagesTest\');

        [rows, cols, cameraPoints, projectorPoints, ~] = decodeMPS(imagesCameraDirectoryPath, params);
        
        tformProjector2Camera = fitgeotrans(projectorPoints, cameraPoints, 'projective');
        homography_cell{i,j} = tformProjector2Camera; 
        
        % Resizing the image so it would fit with the estimated homography
        reconstructionResize = imresize(reconstruction,[params.screenHeight, params.screenWidth]);
        RecImageCameraView = imwarp(reconstructionResize, tformProjector2Camera);

        reconstruction_cell{i,j} = RecImageCameraView;
        
        if debug
            figure,
            colormap(gray);
            imagesc(RecImageCameraView)
            pause;
        end
        
        close all
        disp(['Process for MPS image ',num2str(i),'-',num2str(j),' done.'])
        
        toc
    end 
    
    % skipping check
    if check <= size(cutImage,1)*size(cutImage,2)*5 && g == 2
        continue
    end
    
end
    save('result_test.mat','reconstruction_cell','RecImageCameraView','homography_cell','image_cell');
    
% Loading the result if homography is already calculated
else
    load('result_test.mat','reconstruction_cell','RecImageCameraView','homography_cell','image_cell');
end
