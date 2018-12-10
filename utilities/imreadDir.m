function [ images ] = imreadDir( directoryPath, fileExtension, varargin )
%imreadDir Load all images of specified filetype from specified directory
%   imreadDir loads all images of specified filetype from specified
%   directory to a cell array
%   Optional arguments specify image cropping and image rotation

p=inputParser;

p.addRequired('directoryPath', @isstr);
p.addRequired('fileExtension', @isstr);
p.addParameter('Rotate', 0, @isnumeric);
p.addParameter('Crop', [], @isnumeric);
p.addParameter('Plot', 0, @isnumeric);

p.parse(directoryPath, fileExtension, varargin{:});

% get all files of type fileExtension located in directoryPath
files = dir(fullfile(directoryPath, strcat('*', fileExtension)));
% images=cell(1, length(files));


% load files
for i=1:length(files)
    fileName = files(i).name;
    
    if(~isempty(p.Results.Crop))
        image = im2double(imread(strcat(directoryPath, fileName), 'PixelRegion', {[p.Results.Crop(2),p.Results.Crop(2)+p.Results.Crop(4)], [p.Results.Crop(1),p.Results.Crop(1)+p.Results.Crop(3)]}));
%         image=im2double(imread(strcat(directoryPath, fileName)));
%         image=imcrop(image, p.Results.Crop);
    else
        image=im2double(imread(strcat(directoryPath, fileName)));
    end
    
    if(size(image,3) > 1)
        image = rgb2gray(image);
    end
    
    if(p.Results.Rotate~=0)
        image=imrotate(image, p.Results.Rotate);
    end
    
    if(p.Results.Plot)
        figure(1)
        imagesc(image)
        drawnow
    end
    
    images(:,:,i)=(image);
    disp(['File No: ' num2str(i, '%02d')]);
end

end

