%% Jeffrey Wong | ECE-310 | HW #5

clear
close all
clc

%% Problem 2- Discrete Approximations of the Laplacian

% Part a

% See attached file for derivation

% Sobel operators
hx = [1, 0,-1 ; 2, 0, -2; 1, 0, -1]/4;
hy = [-1, -2, -1; 0, 0, 0; 1, 2, 1]/4;

% Laplacians

hLap = [1, 4, 1; 4, -20, 4; 1, 4, 1]/6;
hLapSob = conv2(hx,hx)+conv2(hy,hy);

% The result is zero-phase due to the values being purely real and being
% symmetric within each row and column

% Part b

[HLap,f1,f2] = freqz2(hLap);
[HLapSob,f1Sob,f2Sob] = freqz2(hLapSob);

% HLap and HLapSob are both nonpositive for all values

% Part c

figure
contour(f1,f2,HLap)
title('Contour Plot of HLap');
xlabel('Frequency in x')
ylabel('Frequency in y')

figure
surf(f1,f2,HLap)
title('Surface Plot of HLap');
xlabel('Frequency in x')
ylabel('Frequency in y')

figure
contour(f1Sob,f2Sob,HLapSob)
title('Contour Plot of HLapSob');
xlabel('Frequency in x')
ylabel('Frequency in y')

figure
surf(f1Sob,f2Sob,HLapSob)
title('Surface Plot of HLapSob');
xlabel('Frequency in x')
ylabel('Frequency in y')

% HLap has larger values near the edges and zero at the center, acting as 
% a highpass filter, but HLapSob acts as a bandpass filter instead
% Both filters appear to be relatively isotropic

% Part d

load("LilyImg.mat")
load("Rodanimg.mat")

figure
image(Lilyx)
title('Grayscale image of Lily');
colormap("gray")

figure
image(Rodanx)
title('Grayscale image of Rodan');
colormap("gray")

LilyLap = filter2(hLap, Lilyx);
LilyLapSob = filter2(hLapSob, Lilyx);
RodanLap = filter2(hLap, Rodanx);
RodanLapSob = filter2(hLapSob, Rodanx);

figure
image(LilyLap)
title('Laplacian applied to Grayscale image of Lily');
colormap("gray")

figure
image(RodanLap)
title('Laplacian applied to Grayscale image of Rodan');
colormap("gray")

figure
image(LilyLapSob)
title('Sobel approximate Laplacian applied to Grayscale image of Lily');
colormap("gray")

figure
image(RodanLap)
title('Sobel approximate Laplacian applied to Grayscale image of Rodan');
colormap("gray")

%% Problem 3- Upsampling in 2D

% Part a- See below

% Part b

test1 = [1,2,3;4,5,6];
test2 = upsamp2(test1);
LilyUpsamp = upsamp2(Lilyx);
RodanUpsamp = upsamp2(Rodanx);

figure
image(LilyUpsamp)
title('Upsampled Grayscale image of Lily');
colormap("gray")

figure
image(RodanUpsamp)
title('Upsampled Grayscale image of Rodan');
colormap("gray")

% Part c

LilyFFT = abs(fftshift(fft2(Lilyx)));
RodanFFT = abs(fftshift(fft2(Rodanx)));

figure
image(LilyFFT)
title('Magnitude of 2D DFT of Lily');
colormap("gray")

figure
image(RodanFFT)
title('Magnitude of 2D DFT of Rodan');
colormap("gray")

LilyUpsampFFT = abs(fftshift(fft2(LilyUpsamp)));
RodanUpsampFFT = abs(fftshift(fft2(RodanUpsamp)));

figure
image(LilyUpsampFFT)
title('Magnitude of 2D DFT of Upsampled Lily');
colormap("gray")

figure
image(RodanUpsampFFT)
title('Magnitude of 2D DFT of Upsampled Rodan');
colormap("gray")

% The upsampled spectrum features imaging distortion

% Part a

function result = upsamp2(A)
   [lenX, lenY] = size(A);
   result = zeros(2*lenX, 2*lenY);
   result(2:2:end, 2:2:end) = A;
end
