debug = 1;

if debug
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
        figure(1122)
        imshow(img);
        hold on;
        rectangle('Position',...
            [roi_wide(1),roi_wide(2),roi_wide(3),roi_wide(4)],...
            'EdgeColor','r','LineWidth',2 )
        drawnow
        hold off;
    end

disp('Found ROI...')
end

roi_wide = round(roi_wide);

%% Import a tiled image and cut it into single pixel-view images.
close all

TiledImagePath = '.\data\Scenes\WhiteScene_all\blackImagesTiled\';
SaveDirPath = '.\data\Scenes\WhiteScene_all\blackImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all images
for g = 1:size(imageFilenames, 2)
% for g = 3:3
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    % For all rows
    % for i = 1:15
    for i = 7:7
        % For all columns
        % for j= 1:15
        for j = 7:7
           % cell_cut{g,i,j} = CutTiled_FixMPS_f(image, roi_wide, i, j, 0);
           A = CutTiled_f(image, roi_wide, i, j, 0);
           % for testing purposes
           imwrite(A, [SaveDirPath, imageFilenames{1,g}]);
        end
    end     
end

disp('Done!')

%%
TiledImagePath = '.\data\Scenes\WhiteScene_all\whiteImagesTiled\';
SaveDirPath = '.\data\Scenes\WhiteScene_all\whiteImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all images
for g = 1:size(imageFilenames, 2)
%for g = 2:2
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    % For all rows
    % for i = 1:15
    for i = 7:7
        % For all columns
        % for j= 1:15
        for j = 7:7
           % cell_cut{g,i,j} = CutTiled_FixMPS_f(image, roi_wide, i, j, 0);
           A = CutTiled_f(image, roi_wide, i, j, 0);
           % for testing purposes
           imwrite(A, [SaveDirPath, imageFilenames{1,g}]);
        end
    end     
end

disp('Done!')

%% 
TiledImagePath = '.\data\Scenes\WhiteScene_all\mpsImagesTiled\';
SaveDirPath = '.\data\Scenes\WhiteScene_all\mpsImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all images
for g = 1:size(imageFilenames, 2)
%for g = 2:2
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    % For all rows
    % for i = 1:15
    for i = 7:7
        % For all columns
        % for j= 1:15
        for j = 7:7
           % cell_cut{g,i,j} = CutTiled_FixMPS_f(image, roi_wide, i, j, 0);
           A = CutTiled_f(image, roi_wide, i, j, 0);
           % for testing purposes
           imwrite(A, [SaveDirPath, imageFilenames{1,g}]);
        end
    end     
end

disp('Done!')

%%
TiledImagePath = '.\data\Scenes\FerScene\mpsImagesTiled\';
SaveDirPath = '.\data\Scenes\FerScene\mpsImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all images
for g = 1:size(imageFilenames, 2)
%for g = 2:2
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    % For all rows
    % for i = 1:15
    for i = 7:7
        % For all columns
        % for j= 1:15
        for j = 7:7
           % cell_cut{g,i,j} = CutTiled_FixMPS_f(image, roi_wide, i, j, 0);
           A = CutTiled_f(image, roi_wide, i, j, 0);
           % for testing purposes
           imwrite(A, [SaveDirPath, imageFilenames{1,g}]);
        end
    end     
end

disp('Done!')
 
%%
TiledImagePath = '.\data\Scenes\FerScene\hadamardImagesTiled\';
SaveDirPath = '.\data\Scenes\FerScene\hadamardImagesTest\';

dirList = dir(TiledImagePath);
isFile = ~[dirList.isdir];
imageFilenames = {dirList(isFile).name};

% For all images
for g = 1:size(imageFilenames, 2)
%for g = 2:2
    image = imread([TiledImagePath, imageFilenames{1,g}]);
    % For all rows
    % for i = 1:15
    for i = 7:7
        % For all columns
        % for j= 1:15
        for j = 7:7
           % cell_cut{g,i,j} = CutTiled_FixMPS_f(image, roi_wide, i, j, 0);
           A = CutTiled_f(image, roi_wide, i, j, 0);
           % for testing purposes
           imwrite(A, [SaveDirPath, imageFilenames{1,g}]);
        end
    end     
end

disp('Done!')

%% For homography testing
% TiledImagePath = '.\data\Scenes\FerScene\whiteImagesTiled\';
% SaveDirPath = '.\data\Scenes\FerScene\whiteImagesTest\';
% 
% dirList = dir(TiledImagePath);
% isFile = ~[dirList.isdir];
% imageFilenames = {dirList(isFile).name};
% 
% cell2 = cell(15,15);
% 
% % For all images
% % for g = 1:size(imageFilenames, 2)
% for g = 2:2
%     image = imread([TiledImagePath, imageFilenames{1,g}]);
%     % For all rows
%     % for i = 1:15
%     for i = 1:15
%         % For all columns
%         % for j= 1:15
%         for j = 1:15
%            % cell_cut{g,i,j} = CutTiled_FixMPS_f(image, roi_wide, i, j, 0);
%            A = CutTiled_f(image, roi_wide, i, j, 0);
%            cell2{i,j} = A;
%            % for testing purposes
% %            imwrite(A, [SaveDirPath, imageFilenames{1,g}]);
%         figure(1)
%         imagesc(cell2{i,j});
%         axis image
%         waitforbuttonpress
%         end
%     end     
% end
% 
% disp('Done!')


% ========================================================================


% %% Saving test images into a folder, this will not be done in the pipeline
% 
% 
% for g = 1:size(imageFilenames, 2)
% %for g = 1:1
%     % For all rows
%     % for i = 1:15
%     for i = 7:7
%         % For all columns
%         % for j= 1:15
%         for j = 7:7
%            imwrite((cell_cut{i,j}), [SaveDirPath, imageFilenames{1,g}]);
%         end
%     end     
% end
% 
% disp('Done!')
% 
% %%
% debug = 1;
% if debug
%     figure()
%     imshow(image);
%     title('Tiled image');
% end
% % 
% % % i represents row in tiled image
% % i = 4;
% % % j represents column in tiled image
% % j = 12;
% 
% % adding for loops for testing
% for i = 1:15
%     for j = 1:15
%         
% % ROI data should be calculated before in the script or provided to the
% % function in a form of an argument. The cut area needs to be calculated
% % with the indexes of the tiled image i and j
% roi_rect(1) = 1 + roi_wide(3)*(j-1);
% roi_rect(2) = 1 + roi_wide(4)*(i-1);
% roi_rect(3) = roi_wide(3);
% roi_rect(4) = roi_wide(4);
% 
% [croppedImage,rect2] = imcrop(image, roi_rect);
% 
% if debug
% %     figure()
% %     imshow(croppedImage);
% %     title('Cropped image');
% 
%     figure(5)
%     imshow(image);
%     hold on;
%     rectangle('Position', [rect2(1),rect2(2),rect2(3),rect2(4)],'EdgeColor','r','LineWidth',2 );
%     title('Showing a position that was cropped')
%     drawnow
% 
% end
%     end
% end
% 
% 
