function [imagesCamera] = captureImages(imagesProjectorDirectoryPath, imagesCameraDirectoryPath, params)
%capturePhaseShift Summary of this function goes here
%   Detailed explanation goes here


dirList = dir(imagesProjectorDirectoryPath);
isFile = ~[dirList.isdir];
imageFilenames = sort({dirList(isFile).name});

fullscreen(zeros([params.screenHeight, params.screenWidth]), params.screenID)

if(exist(imagesCameraDirectoryPath,'dir')~=7)
    mkdir(imagesCameraDirectoryPath)
end

for i = 1:params.captureEveryNthImage:size(imageFilenames, 2)
    display(['Capturing:     (',num2str(i),'/',num2str(size(imageFilenames, 2)),') ', imageFilenames{i}])
    
    tic
    
    imageProjector = imread([imagesProjectorDirectoryPath, imageFilenames{i}]);
    
    if(i==1)
        pause(1.3*params.delayTimeSlideShow)
    else
        pause(params.delayTimeSlideShow)
    end
    
    fullscreen(imageProjector, params.screenID);
    
    %     beep
    
    if(i==1)
        pause(1.3*params.delayTimeSlideShow)
    else
        pause(params.delayTimeSlideShow)
    end
    
    [imagesCamera(:,:,i), timeStamp(:,:,i)] = snapshot(params.cam);
    
    
    if(params.writeFlag && strcmp(params.writeMode, 'seq'))
        
        [~,name,~] = fileparts(imageFilenames{i});
        display(['Writing:     (',num2str(i),'/', num2str(size(imageFilenames, 2)),') ', imagesCameraDirectoryPath, name, '.tiff'])
        
        if(params.cropFlag)
            image = imcrop(imagesCamera(:,:,i), round([params.roi(1), 0, params.roi(3), size(imagesCamera(:,:,i), 1)]));
        end
        
        imwrite(image, [imagesCameraDirectoryPath, name, '.tiff']);
    end
    
    if(i==1)
        pause(1.3*params.delayTimeCamera)
    else
        pause(params.delayTimeCamera)
    end
    
    toc
    display(newline)
    
end

closescreen

if(params.writeFlag && strcmp(params.writeMode, 'batch'))
    for i = 1:size(imageFilenames, 2)
        [~,name,~] = fileparts(imageFilenames{i});
        display(['Writing:     (',num2str(i),'/', num2str(size(imageFilenames, 2)),') ', imagesCameraDirectoryPath, name, '.tiff'])
        
        if(params.cropFlag)
            image = imcrop(imagesCamera(:,:,i), round([params.roi(1), 0, params.roi(3), size(imagesCamera(:,:,i), 1)]));
        end
        
        imwrite(image, [imagesCameraDirectoryPath, name, '.tiff']);
        
        
    end
end


end

