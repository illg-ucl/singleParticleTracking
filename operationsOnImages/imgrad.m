function g = imgrad(image)
%
% Created by Isabel Llorente-Garcia, August 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% calculate the gradient of an image (single frame).
h=fspecial('sobel'); % filter mask to calculate gradient: h=[1 2 1; 0 0 0; -1 -2 -1] and h'=[1 0 -1; 2 0 -2; 1 0 -1].
% Calculate gradiente image:
grady = imfilter(image,h,'replicate'); % vertical gradient.
gradx = imfilter(image,h','replicate'); % horizontal gradient.
g = sqrt(grady.^2+gradx.^2); % gradient image.