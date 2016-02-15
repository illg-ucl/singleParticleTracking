function y  = normalisedGaussian(x,mean,sigma)
%
% Created by Isabel Llorente-Garcia, June 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Normalised Gaussian function.
%
% Input x can be a value or a vector.
% Inputs mean and sigma are values.
%
% Example:
% x = (0:0.01:20);
% y = normalisedGaussian(x,4,1);
% plot(x,y)
%
y = 1/(sqrt(2*pi)*sigma)*exp(-(x-mean).^2/(2*sigma^2));
