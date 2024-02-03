%% Jeffrey Wong | ECE-310 | HW #3

clear
close all
clc

%% Problem 5- Wavelets and Filter Banks

% Part a
[h0_1,h1_1,f0_1,f1_1] = wfilters('db1');

w = linspace(0,pi,1e4);
h0_1_resp = freqz(h0_1,1,w);
h1_1_resp = freqz(h1_1,1,w);

figure
hold on
legend
plot(w,abs(h0_1_resp),'DisplayName',"H_0(w)");
plot(w,abs(h1_1_resp),'DisplayName',"H_1(w)");
title('Frequency Response of Haar Filters');
xlabel('Frequency (rad)');
ylabel('Magnitude');

% See attached file for E(z) and verification of paraunitary PR

% Part b
P_haar = abs(h0_1_resp).^2 + abs(h1_1_resp).^2;
P_haar_avg = mean(P_haar); % The average power should be 2
P_haar_range = max(P_haar)-min(P_haar); % Power seems to have minimal variation

% Part c
[h0_5,h1_5,f0_5,f1_5] = wfilters('db5');

% 1) The dbN filters seem to have a length of 2N. 
% Generally, f_0[n] = h_0[2N - n -1], h_1[n] = -h_0[2N - n - 1] for even n and
% h_0[2N - n - 1] for odd n, and f_1[n] = h_0[n] for even n and -h_0[n] for odd n

% 2)

E = zeros(2,2,5); 
% "Even" components have odd indices due to 1-indexing, and vice versa
% Strictly speaking these values of E represent the coefficients of the 
% polyphase components only, relevant for phase 3
e11 = h0_5(1:2:9);
E(1,1,:) = e11;
e12 = h0_5(2:2:10);
E(1,2,:) = e12;
e21 = h1_5(1:2:9);
E(2,1,:) = e21;
e22 = h1_5(2:2:10);
E(2,2,:) = e22;

% 3)

% Note that as per #3, fliplr gives z^-(L-1) times the paraconjugate since
% entries are all entirely real

% This means that polynomial multiplication via conv would give 9 entries- 
% 1 for z^4 all the way down to 9 for z^-4

result = zeros(2,2,9);

result(1,1,:) = conv(fliplr(e11),e11) + conv(fliplr(e21),e21);
result(1,2,:) = conv(fliplr(e11),e12) + conv(fliplr(e21),e22);
result(2,1,:) = conv(fliplr(e12),e11) + conv(fliplr(e22),e21);
result(2,2,:) = conv(fliplr(e12),e12) + conv(fliplr(e22),e22);

% result(:,:,5) corresponds to our constant (z^0) terms- should equal
% identity. All other layers are essentially zero- our result is the
% identity matrix, which means we satisfy PR paraunity property.

% 4) Since ~E * E = I, R(z) = ~E(z), meaning N = 1. M = 2, so end-to-end
% delay is MN - 1 = 2 - 1 = 1

% 5)

h0_5_resp_1 = freqz(h0_5,1,w);
h1_5_resp_1 = freqz(h1_5,1,w);

figure
hold on
legend
plot(w,abs(h0_5_resp_1),'DisplayName',"H_0(w)");
plot(w,abs(h1_5_resp_1),'DisplayName',"H_1(w)");
xlabel('Frequency (rad)');
ylabel('Magnitude');
title('Frequency Response of Order 5 Daubechies Filters');

% 6) Isn't this the same as 3)?

% 7)

h0_5_resp_2 = freqz(h0_5,1,2*w);
h1_5_resp_2 = freqz(h1_5,1,2*w);
h0_5_resp_4 = freqz(h0_5,1,4*w);
h1_5_resp_4 = freqz(h1_5,1,4*w);

g_0 = h0_5_resp_1.*h0_5_resp_2.*h0_5_resp_4;
g_1 = h0_5_resp_1.*h0_5_resp_2.*h1_5_resp_4;
g_2 = h0_5_resp_1.*h1_5_resp_2;
g_3 = h1_5_resp_1;

figure
hold on
legend
plot(w,abs(g_0),'DisplayName',"G_0(w)");
plot(w,abs(g_1),'DisplayName',"G_1(w)");
plot(w,abs(g_2),'DisplayName',"G_2(w)");
plot(w,abs(g_3),'DisplayName',"G_3(w)");
xlabel('Frequency (rad)');
ylabel('Magnitude');
title('Frequency Response of Tree Structure Channels');

% Passband widths decrease as frequency

% 8)

P_5 = abs(g_0).^2 / 8 + abs(g_1).^2 / 8 + abs(g_2).^2 / 4 + abs(g_3).^2 / 2;
P_5_avg = mean(P_5); % The average power should be 1
P_5_range = max(P_5)-min(P_5); % Power seems to have minimal variation (Less than 10^-10)

% Since the filter bank is maximally decimated, sum(1/Mk) is 1.