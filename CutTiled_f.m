function [croppedImage] = CutTiled_f(image, roi_wide, rowNo, colNo, debug)

% Import a tiled image and cut it into single pixel-view images.
        
% ROI data should be calculated before in the script or provided to the
% function in a form of an argument. The cut area needs to be calculated
% with the indexes of the tiled image i and j
roi_rect(1) = 1 + roi_wide(3)*(colNo-1);
roi_rect(2) = 1 + roi_wide(4)*(rowNo-1);
roi_rect(3) = roi_wide(3);
roi_rect(4) = roi_wide(4);

[croppedImage,rect2] = imcrop(image, roi_rect);

if debug

    figure(5)
    imshow(image);
    hold on;
    rectangle('Position', [rect2(1),rect2(2),rect2(3),rect2(4)],'EdgeColor','r','LineWidth',2 );
    title('Showing a position that was cropped')
    drawnow

end

end

