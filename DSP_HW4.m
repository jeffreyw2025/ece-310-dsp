%% Jeffrey Wong | ECE-310 | HW #4

clear
close all
clc

%% Problem 4- SOS Representation of Filters

% Part a

[z,p,k] = ellip(4,1.5,30,[0.3 0.6], "bandpass"); % Prototype order 4 gives us bandpass order 8
[b,a]= zp2tf(z,p,k);

w = linspace(0,pi,1e4);
h = freqz(b,a,w);

figure
plot(w,20*log10(abs(h)));
xlabel('Frequency (rad)');
ylabel('Gain (db)');
title('Magnitude Response of Bandpass Filter');
axis([0 pi -50 2]);

% Part b

[sosUp, gUp] = zp2sos(z,p,k,"up","inf");
[sosDown, gDown] = zp2sos(z,p,k,"down","inf");

% Part c

[zu1,pu1,ku1] = tf2zpk(sosUp(1,1:3),sosUp(1,4:6)); 
[zu2,pu2,ku2] = tf2zpk(sosUp(2,1:3),sosUp(2,4:6)); 
[zu3,pu3,ku3] = tf2zpk(sosUp(3,1:3),sosUp(3,4:6)); 
[zu4,pu4,ku4] = tf2zpk(sosUp(4,1:3),sosUp(4,4:6));

zUp = [abs(zu1) abs(zu2) abs(zu3) abs(zu4)];
pUp = [abs(pu1) abs(pu2) abs(pu3) abs(pu4)];

[zd1,pd1,kd1] = tf2zpk(sosDown(1,1:3),sosDown(1,4:6)); 
[zd2,pd2,kd2] = tf2zpk(sosDown(2,1:3),sosDown(2,4:6)); 
[zd3,pd3,kd3] = tf2zpk(sosDown(3,1:3),sosDown(3,4:6)); 
[zd4,pd4,kd4] = tf2zpk(sosDown(4,1:3),sosDown(4,4:6));

zDown = [abs(zd1) abs(zd2) abs(zd3) abs(zd4)];
pDown = [abs(pd1) abs(pd2) abs(pd3) abs(pd4)];

% Poles are ordered in magnitude correctly, but zeroes are all 1?

% Part d

aUp1 = sosUp(1,4:6);
hUp1 = freqz(gUp,aUp1,w);
bUp1 = gUp*sosUp(1,1:3);
aUp2 = conv(aUp1, sosUp(2,4:6));
hUp2 = freqz(bUp1,aUp2,w);
bUp2 = conv(bUp1, sosUp(2,1:3));
aUp3 = conv(aUp2, sosUp(3,4:6));
hUp3 = freqz(bUp2,aUp3,w);
bUp3 = conv(bUp2, sosUp(3,1:3));
aUp4 = conv(aUp3, sosUp(4,4:6));
hUp4 = freqz(bUp3,aUp4,w);

figure
hold on
legend
plot(w,20*log10(abs(hUp1)),'DisplayName',"F_1(w)");
plot(w,20*log10(abs(hUp2)),'DisplayName',"F_2(w)");
plot(w,20*log10(abs(hUp3)),'DisplayName',"F_3(w)");
plot(w,20*log10(abs(hUp4)),'DisplayName',"F_4(w)");
xlabel('Frequency (rad)');
ylabel('Gain (dB)');
title('Magnitude Response of F filters for up ordering');
axis([0 pi -50 2]);

aDown1 = sosDown(1,4:6);
hDown1 = freqz(gDown,aDown1,w);
bDown1 = gDown*sosDown(1,1:3);
aDown2 = conv(aDown1, sosDown(2,4:6));
hDown2 = freqz(bDown1,aDown2,w);
bDown2 = conv(bDown1, sosDown(2,1:3));
aDown3 = conv(aDown2, sosDown(3,4:6));
hDown3 = freqz(bDown2,aDown3,w);
bDown3 = conv(bDown2, sosDown(3,1:3));
aDown4 = conv(aDown3, sosDown(4,4:6));
hDown4 = freqz(bDown3,aDown4,w);

figure
hold on
legend
plot(w,20*log10(abs(hDown1)),'DisplayName',"F_1(w)");
plot(w,20*log10(abs(hDown2)),'DisplayName',"F_2(w)");
plot(w,20*log10(abs(hDown3)),'DisplayName',"F_3(w)");
plot(w,20*log10(abs(hDown4)),'DisplayName',"F_4(w)");
xlabel('Frequency (rad)');
ylabel('Gain (dB)');
title('Magnitude Response of F filters for down ordering');
axis([0 pi -50 2]);

% Part e

% The ordering of the denominator polynomials was verified to be reversed
% by inspection. 

numRatio1 = sosUp(1,1:3)./sosDown(4,1:3);
numRatio2 = sosUp(2,1:3)./sosDown(3,1:3);
numRatio3 = sosUp(3,1:3)./sosDown(2,1:3);
numRatio4 = sosUp(4,1:3)./sosDown(1,1:3);

% These should be within floating-point error of 0
varNumRatio1 = abs(max(numRatio1) - min(numRatio1));
varNumRatio2 = abs(max(numRatio2) - min(numRatio2));
varNumRatio3 = abs(max(numRatio3) - min(numRatio3));
varNumRatio4 = abs(max(numRatio4) - min(numRatio4));

%% Problem 5- Sensitivity Properties of Parallel Allpass Realizations

% Part a

% Original filter coefficents
b0 = [0.1336 0.0563 0.0563 0.1336];
a0 = [1 -1.5055 1.2630 -0.3778];

ybq0 = fi(b0, 1, 5, 3); % We need 1 integer bit to represent the 1.5 value and 1 sign bit, can only use 3 fractional bits
bq0 = ybq0.data; % Quantized numerator of original tf

yaq0 = fi(a0, 1, 5, 3);
aq0 = yaq0.data; % Quantized denominator of original tf

% Allpass filter coefficents
bA1 = [-0.4954 1];
aA1 = [1 -0.4954];
bA2 = [0.7626 -1.0101 1];
aA2 = [1 -1.0101 0.7626];

ybqA1 = fi(bA1, 1, 5, 3);
bqA1 = ybqA1.data;
yaqA1 = fi(aA1, 1, 5, 3);
aqA1 = yaqA1.data;
ybqA2 = fi(bA2, 1, 5, 3);
bqA2 = ybqA2.data;
yaqA2 = fi(aA2, 1, 5, 3);
aqA2 = yaqA2.data;

% Part b

% See attached file for computations

sumb0 = sum(b0); % For some reason these sums do not match so plugging in z = 1 does not give exactly 1 for the gain
suma0 = sum(a0);

% Error only happens in the quantized origninal, where H_Q0(1) = 2/3

errQ0 = abs(20*log10(2/3)); % log10(1) = 0 so we ignore it

% Part c

h0 = freqz(b0, a0, w);
magh0 = 20*log10(abs(h0));
hq0 = freqz(bq0, aq0, w);
maghq0 = 20*log10(abs(hq0));
hdiff0 = abs(h0-hq0);
hqA1 = freqz(0.5*bqA1, aqA1, w); % Allpass factors have a constant 1/2
hqA2 = freqz(0.5*bqA2, aqA2, w); 
hqA = hqA1 + hqA2; % Frequency response is equal to the sum of the frequency responses of the allpass factors
maghqA = 20*log10(abs(hqA));
hdiffA = abs(h0-hqA);

figure
hold on
legend
plot(w,magh0,'DisplayName',"H(w) [Original]");
plot(w,maghq0,'DisplayName',"H_{Q0}(w) [Quantized Orignal]");
plot(w,maghqA,'DisplayName',"H_{QA}(w) [Quantized Allpass]");
xlabel('Frequency (rad)');
ylabel('Gain (dB)');
title('Magnitude Response of Original and Quantized Filters');
axis([0 pi -40 2]);

% Part d

% Max differences between ideal and quantized original and allpass filters
maxhdiff0 = max(abs(hdiff0));
maxhdiffA = max(abs(hdiffA)); % Deviation should be lower for the allpass filter

% The quantized original filter does not meaningfully "ripple" at all, 
% while the quantized allpass filter closely follows the ideal filter's
% equiripple behavior in the passband

% Compared to the ideal filter, the quantized original filter's stopband 
% only has ~15 dB attenuation in the stopband, while the quantized allpass
% filter's stopband actually achieves more attenuation than the original at
% about 21.4 dB

% Part e

% Figure out how to compute group delay properly (mostly what units)

% Phase responses of filters (in radians)
phase0 = unwrap(angle(h0(1:3001)));
phaseq0 = unwrap(angle(hq0(1:3001)));
phaseqA = unwrap(angle(hqA(1:3001)));


grp0 = -(diff(phase0)/(w(2)-w(1))); % Passband is 3/10 normalized to Nyq. Bandwith, so we take the first 3/10 of samples
grpq0 = -(diff(phaseq0)/(w(2)-w(1)));
grpqA = -(diff(phaseqA)/(w(2)-w(1)));

figure
hold on
legend
plot(w(1:3000),grp0,'DisplayName',"H(w) [Original]");
plot(w(1:3000),grpq0,'DisplayName',"H_{Q0}(w) [Quantized Orignal]");
plot(w(1:3000),grpqA,'DisplayName',"H_{QA}(w) [Quantized Allpass]");
xlabel('Frequency (rad)');
ylabel('Group Delay');
title('Group Delay of Original and Quantized Filters in Passband');

grpdevq0 = max(abs(grp0-grpq0));

grpdevqA = max(abs(grp0-grpqA));

% The parallel allpass realization shows much better sensitivity when
% quantized, as it is much closer to the group delay of the non-quantized 
% filter than the quantized original filter.

% Part f

[z0,p0,k0] = tf2zpk(b0, a0);
[zq0,pq0,kq0] = tf2zpk(bq0, aq0);

bqA = 0.5*(conv(bqA1,aqA2) + conv(aqA1,bqA2));
aqA = conv(aqA1,aqA2);
[zqA,pqA,kqA] = tf2zpk(bqA, aqA);


figure
zplane(z0,p0)
title('Pole-Zero Plot of Orignal Filter');

figure
zplane(zq0,pq0)
title('Pole-Zero Plot of Quantized Orignal Filter');

figure
zplane(zqA,pqA)
title('Pole-Zero Plot of Quantized Allpass Filter');

% Note: it seems that the zero at z = -1 does not move for any of the
% quantized filters. The other zeroes did move, with the real part of the
% zeroes increasing in the quantized original filter and decreasing
% slightly in the quantized allpass filter


