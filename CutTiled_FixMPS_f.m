function [croppedImage] = CutTiled_FixMPS_f(image, roi_wide, rowNo, colNo, debug)

% Import a tiled image and cut it into single pixel-view images. Here;
% function is cutting the area 5% smaller than the original ROI in the
% scene. That is beacuse (probably) noise is causing MPS function to
% produce wrong results that make whole devignetting and denoising steps
% not usefull.

% ROI data should be calculated before in the script or provided to the
% function in a form of an argument. The cut area needs to be calculated
% with the indexes of the tiled image i and j
roi_rect(1) = 1 + roi_wide(3)*(colNo-1);
roi_rect(2) = 1 + roi_wide(4)*(rowNo-1);
roi_rect(3) = roi_wide(3);
roi_rect(4) = roi_wide(4);

[croppedImage1,~] = imcrop(image, roi_rect);

% New ROI coordinates are calculated from the cropped image. The process is
% the same, but the starting point is (1,1) and width and length of the
% image are used to calculate new width and length.
% 10% of the width (x dim)
roi_scale(1) = size(croppedImage1,2)*0.1;
% 10% of the length (y dim)
roi_scale(2) = size(croppedImage1,1)*0.1;

roi_MPS(1) = round(1 + roi_scale(1));
roi_MPS(2) = round(1 + roi_scale(2));
roi_MPS(3) = round(size(croppedImage1,2) - 2*roi_scale(1));
% y dim is cut for a little bit more (25%) because of the scene structure
roi_MPS(4) = round(size(croppedImage1,1) - 3*roi_scale(2));

[croppedImage,rect2] = imcrop(croppedImage1, roi_MPS);

if debug
figure(3)
imshow(croppedImage1);
hold on;
rectangle('Position', [rect2(1),rect2(2),rect2(3),rect2(4)],'EdgeColor','r','LineWidth',2 );
hold off;
end

if debug
    
    figure(5)
    imshow(image);
    hold on;
    rectangle('Position', [rect2(1),rect2(2),rect2(3),rect2(4)],'EdgeColor','r','LineWidth',2 );
    title('Showing a position that was cropped')
    drawnow

end

end

