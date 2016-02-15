function g = dftcorr(f,w)
%
% Created by Isabel Llorente-Garcia, August 2011. 
%
% See pg 491, book: "Digital Image Processing using MATLAB". Rafael C.
% Gonzalez, Richard E. Woods and Steven L. Eddins.
%
% dftcorr: 2D correlation in the frequency domain.
% Performs the correlation of a mask, w, with image f. The output g is the
% correlation image, of class double. The output is of the same size as f.
% When, as is generally true in practice, the mask image is much smaller
% than f, wraparound error is negligible if w is padded to size(f).
%


[M, N] = size(f);
f = fft2(f); % Fourier Transform of the image.
w = conj(fft2(w,M,N)); % Conjugate of Fourier Transform of mask image.
g = real(ifft2(w.*f)); % Inverse Fourier transform of the product of the last two.
