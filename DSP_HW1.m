%% Jeffrey Wong | ECE-310 | HW #1

clear
close all
clc

%% Problem 2

% See attached file for Problem 1 and Parts a and b for Problem 2.

% Part c
N = 512; % We are performing a 512-pt DFT
M = 500; % Number of samples
fs = 2e7; % Sampling rate
fbin = fs/N; % Bin frequency
fwave = 6e6; % Signal frequency

% Generating the sinewave
k = 0:M-1;
t = 1/(fs) * k;
x = sin(2*pi*fwave*t); % We sample 500 times at fs then pad with zeroes
xrect = x / sqrt(M); % A rectangular window has value 1 over its support so it has a total energy equal to the length

% Plotting the DFT of the rectangular window
Xrect = fft(xrect,N);
irect = -N/2:N/2-1;
frect = fs/N .* irect;
figure
subplot(2,1,1)
plot(frect,abs(fftshift(Xrect)))
title('Magnitude of Fourier Transform of Sinewave w/ Rectangular Window')
xlabel('Frequency (Hz)')
ylabel('Linear Magnitude')
subplot(2,1,2)
plot(frect,unwrap(angle(fftshift(Xrect))*180/pi))
title('Phase of Fourier Transform of Sinewave w/ Rectangular Window')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')

% Plotting the DFT of the Chebyshev window
cwin = transpose(chebwin(M,30)); % chebwin generates a column vector, we want a row vector
cwin = cwin / sqrt(sum(cwin.^2)); % The total energy of the Chebyshev window is given by the sum of its index values squared
xcheb = x .* cwin;
Xcheb = fft(xcheb,N);
figure
subplot(2,1,1)
plot(frect,abs(fftshift(Xcheb)))
title('Magnitude of Fourier Transform of Sinewave w/ Chebyshev Window')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
subplot(2,1,2)
plot(frect,unwrap(angle(fftshift(Xcheb))*180/pi))
title('Phase of Fourier Transform of Sinewave w/ Chebyshev Window')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')


izoom = 144:164; % Peak value k_0 is expected to be 6/20 * 512 = 153.6 ~ 154
figure
legend
title('Zoomed Plot of Magnitude of Fourier Transform of Sinewave')
xlabel('Bin Index')
ylabel('Magnitude')
hold on
plot(izoom, abs(Xrect(145:165)),'DisplayName',"Rectangular Window") % Note bins start at 0 while vectors are 1-indexed, so adding 1 was necessary
plot(izoom, abs(Xcheb(145:165)),'DisplayName',"Chebyshev Window")
hold off

% Part d

xnoisy = addNoise(x,20);
arect = (xnoisy)/sqrt(M);
Arect = abs(fft(arect,512));
Arect = 20*log10(Arect); % Gets the magnitude of Arect
figure
plot(frect,fftshift(Arect))
title('Fourier Transform of Noisy Sinewave w/ Rectangular Window')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

acheb = (xnoisy).*cwin;
Acheb = abs(fft(acheb,512));
Acheb = 20*log10(Acheb); % Gets the magnitude of Acheb

figure
plot(frect,fftshift(Acheb))
title('Fourier Transform of Noisy Sinewave w/ Chebyshev Window')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

fzoom = izoom .* (fs/N);

figure
hold on
legend
title('Zoomed Plot of Magnitude of FT of Noisy Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
plot(fzoom,Arect(145:165),'DisplayName',"Rectangular Window")
plot(fzoom,Acheb(145:165),'DisplayName',"Chebyshev Window")
hold off

%% Problem 3

% Part a

bpm = 100;
fs = 44100;
notes = 63; % Number of quarter notes to generate- See part b for how this number was calculated
ndB = -40; 
N = 10^(ndB/20)*randn(1,notes); % Noise vector

% We generate the tones included by selecting a tone to be excluded then
% taking the complement
excludedTones = randi(3,1,notes);
[U,V] = meshgrid(excludedTones,[1;2;3]); % Cols of U denote the excluded note, Rows of V denote tone 1, 2, or 3
included = ~(U == V); % U == V gives the excluded tones at each notes

g4 = generateBeats(notes, bpm, fs, 392, included(1,:));
a4 = generateBeats(notes, bpm, fs, 440, included(2,:));
d5 = generateBeats(notes, bpm, fs, 587.33, included(3,:));
combinedBeats = addNoise(g4 + a4 + d5, 40);

% Part b

N_0 = 2^15; % A bin width of 2 Hz requires 22050 bins, 32768 or 2^15 is the next largest power of 2
% The total number of samples is given by N + (M-1)L = N_0 + 99(N_0/2) = 1654784 for our situation
% Since each quarter note has (60/100)*44100 = 26460 samples, we need at least 62.5 notes, we round to 63
[pdg,fwel] = pwelch(combinedBeats,hamming(N_0),N_0/2,N_0,fs);

% Part c

figure
plot(fwel,20*log10(pdg)) % The dB plot shows the tail a bit better
title('Welch Periodogram of Note Samples')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

% Lower and upper bound for zoomed plot
lb = find(fwel<300,1,"last");
ub = find(fwel>600,1,"first");

figure
plot(fwel(lb:ub),20*log10(pdg(lb:ub)))
title('Zoomed-in Welch Periodogram of Note Samples from 300-600 Hz')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-50 -10])

% Part d

% Note from part b that each quarter note has 26460 samples, so each eighth
% note would have 13230 samples

figure
spectrogram(combinedBeats,hamming(N_0),13230,N_0);

% Function for part a

function bVector = generateBeats(beats, bpm, fs, ftone, included)
    beatlength = 60/bpm; % Length of beat in seconds is reciporical of beats per second
    beatsize = round(beatlength * fs); % Number of sampled points per beat
    tsamp = (0:beatsize-1)./fs; % Round to ensure we have an integer, we sample every 1/fs seconds
    basebeat = sin((2*pi*ftone).*tsamp);
    beats = repmat(basebeat,beats,1);
    bMatrix = beats .* repmat(included', 1, beatsize);
    bVector = reshape(bMatrix',1,[]); % Transform the grid into a vector. Note that matrix has to be transposed
    % because reshape processes by columns first
end

function newVector = addNoise(vector, SNR)
    vsize = length(vector);
    noise = rms(vector)*10^(-SNR/20)*randn(1,vsize); % rms gives the standard deviation of the signal
    newVector = vector + noise;
end