clear all
close all
 
% Finding ROI in image

ImagesDirPath = '.\data\imagesCamera\Tiff\';

% variable debug for displaying images
debug = 1;

%% Section used for processing whole folder
% dirList = dir(ImagesDirPath);
% isFile = ~[dirList.isdir];
% imageFilenames = {dirList(isFile).name};

%% Section for processing a single picture

    testImg = imread([ImagesDirPath, 'slika_5.tiff']);
    
    if debug == 1
        figure()
        imshow(testImg);
        title('Original Image');
    
        %%%%
    end
    
    %[counts,x] = imhist(im_cell{i,j},1024); 
    [counts,x] = imhist(testImg,1024); 
    T = otsuthresh(counts); 
    testImgBin = imbinarize(testImg, 0.75*T);
    %im_cell_bin{i,j} = imbinarize(im_cell{i,j}, 0.75*T);
    
    % find area of all contours and corresponding bboxes; define the ROI as the 
    % bbox of the largest area contour 
    s = regionprops(testImgBin, 'area', 'boundingBox'); 
    [~, maxAreaIdx] = max(cell2mat({s.Area}));
                roi = s(maxAreaIdx).BoundingBox; 
    
    if debug == 1
        figure()
        imshow(testImgBin);
        hold on;
        rectangle('Position', [roi(1),roi(2),roi(3),roi(4)],'EdgeColor','r','LineWidth',2 )
        %plot(x, y, 'r.-', 'MarkerSize', 15);
    end