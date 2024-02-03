%% Jeffrey Wong | ECE-310 | HW #2

% Please check command window for certain outputs. Note that program runs
% in filter type order, not part order

clear
close all
clc

%% Problem 2- IIR Filter Design

% See attached file for Problem 1 and Part a of Problem 3

% TODO: Fix axis scaling for analog filter zplanes

% Specifications for all filters
fsamp = 4e7;
fnyq = fsamp/2;
Wp_d = [9 12.5]/20; % Passband frequencies normalized to Nyq. bandwidth
Ws_d = [9.5 12]/20; % Stopband frequencies normalized to Nyq. bandwidth
Wp_a = Wp_d * pi * fsamp; % Passband frequencies as angular analog frequencies
Ws_a = Ws_d * pi * fsamp;
Rp = 1.5; % Passband variation
Rs = 40; % Stopband attenuation
f = linspace(0,fnyq,1000); % Analog and digital frequency vectors for plots 

% Digital Butterworth

% Part a
[n,wn] = buttord(Wp_d,Ws_d,Rp,Rs); % Digital Butterworth Filter Prototype
[z,p,k] = butter(n,wn,"stop");

% Part b
figure
zplane(z,p)
title('Pole-Zero Plot of Digital Butterworth Filter');

% Parts c&d
[b,a] = zp2tf(z,p,k);
h = freqz(b,a,pi*f/fnyq);

displayDigitalFilterInfo(length(p), f, h, "Digital Butterworth") % length(p) gives the degree of the denominator and thus the actual filter

% Analog Butterworth Filter

% Part a
[n,wn] = buttord(Wp_a,Ws_a,Rp,Rs,'s'); % Analog Butterworth Filter Prototype- Note frequencies must be in rad/sec
[z,p,k] = butter(n,wn,"stop",'s');

% Part b

figure
zplane(z,p)
title('Pole-Zero Plot of Analog Butterworth Filter');

% Parts c&d

[b,a] = zp2tf(z,p,k);
h = freqs(b,a,2*pi*f);
displayDigitalFilterInfo(length(p), f, h, "Analog Butterworth")

% Digital Chebyshev I

% Part a
[n,wp] = cheb1ord(Wp_d,Ws_d,Rp,Rs); % Digital Cheb. I Filter Prototype
[z,p,k] = cheby1(n,Rp,wp,"stop");

% Part b
figure
zplane(z,p)
title('Pole-Zero Plot of Digital Chebyshev I Filter');

% Parts c&d
[b,a] = zp2tf(z,p,k);
h = freqz(b,a,pi*f/fnyq);
displayDigitalFilterInfo(length(p), f, h, "Digital Chebyshev I")

% Analog Chebyshev I

% Part a
[n,wp] = cheb1ord(Wp_a,Ws_a,Rp,Rs,'s'); % Analog Butterworth Filter Prototype- Note frequencies must be in rad/sec
[z,p,k] = cheby1(n,Rp,wp,"stop",'s');

% Part b

figure
zplane(z,p)
title('Pole-Zero Plot of Analog Chebyshev I Filter');

% Parts c&d

[b,a] = zp2tf(z,p,k);
h = freqs(b,a,2*pi*f);
displayDigitalFilterInfo(length(p), f, h, "Analog Chebyshev I")

% Digital Chebyshev II

% Part a
[n,ws] = cheb2ord(Wp_d,Ws_d,Rp,Rs); % Digital Cheb. II Filter Prototype
[z,p,k] = cheby2(n,Rs,ws,"stop");

% Part b
figure
zplane(z,p)
title('Pole-Zero Plot of Digital Chebyshev II Filter');

% Parts c&d
[b,a] = zp2tf(z,p,k);
h = freqz(b,a,pi*f/fnyq);
displayDigitalFilterInfo(length(p), f, h, "Digital Chebyshev II")

% Analog Chebyshev II

% Part a
[n,ws] = cheb2ord(Wp_a,Ws_a,Rp,Rs,'s'); % Analog Cheb. II Filter Prototype
[z,p,k] = cheby2(n,Rs,ws,"stop",'s');

% Part b
figure
zplane(z,p)
title('Pole-Zero Plot of Analog Chebyshev II Filter');

% Parts c&d
[b,a] = zp2tf(z,p,k);
h = freqs(b,a,2*pi*f);
displayDigitalFilterInfo(length(p), f, h, "Analog Chebyshev II")

% Digital Elliptic Filter

% Part a
[n,wn] = ellipord(Wp_d,Ws_d,Rp,Rs); % Elliptic Filter Prototype
[z,p,k] = ellip(n,Rp,Rs,wn,"stop");

% Part b
figure
zplane(z,p)
title('Pole-Zero Plot of Digital Elliptic Filter');

% Parts c&d
[b,a] = zp2tf(z,p,k);
h = freqz(b,a,pi*f/fnyq);
displayDigitalFilterInfo(length(p), f, h, "Digital Elliptic")

% Analog Elliptic Filter

% Part a
[n,wn] = ellipord(Wp_a,Ws_a,Rp,Rs,'s'); % Elliptic Filter Prototype
[z,p,k] = ellip(n,Rp,Rs,wn,"stop",'s');

% Part b
figure
zplane(z,p)
title('Pole-Zero Plot of Analog Elliptic Filter');

% Parts c&d
[b,a] = zp2tf(z,p,k);
h = freqs(b,a,2*pi*f);
displayDigitalFilterInfo(length(p), f, h, "Analog Elliptic")

%% Problem 3- FIR Filter Design

% Part a- See attached work for full derivation

delstop = 1e-2;
delpass = (10^(3/40)-1)/(10^(3/40)+1);
dev = [delpass delstop delpass];

% General Variables

fedges = [9e6 9.5e6 12e6 12.5e6];
a = [1 0 1]; % Amplitude in passband and stopband
lpassedge = find((f < 9e6),1,"last");
lstopedge = find((f > 9.5e6),1,"first");
rstopedge = find((f < 12e6),1,"last");
rpassedge = find((f > 12.5e6),1,"first");

% Kaiser Filter

% Part b

[n,Wn,beta,ftype] = kaiserord(fedges,a,dev,fsamp);
window = kaiser(n+1,beta);
b_kaiser = fir1(n,Wn,ftype,window);
disp("Order of Kaiser FIR Filter: " + n);

% Stem plots of Filter coefficients
figure
stem(0:length(b_kaiser)-1,b_kaiser,'ro')
title("Kaiser FIR Filter Coefficients")

figure
zplane(b_kaiser,1)
title('Pole-Zero Plot of Kaiser FIR Filter');
    
% Frequency Response
h = freqz(b_kaiser,1,f,fsamp);
mag = 20*log10(abs(h));
pha = unwrap(angle(h))*180/pi;
    
figure
subplot(2,1,1)
plot(f/1e6,mag);
xlabel('Frequency (MHz)');
ylabel('Gain (db)');
title('Magnitude Response of Kaiser FIR Filter');
grid on;
axis([0 20 -50 2]);
    
subplot(2,1,2)
plot(f/1e6,pha);
xlabel('Frequency (MHz)');
ylabel('Phase (deg)');
title('Phase Response of Kaiser FIR Filter');
grid on;

% Part d

m_passband = [mag(1:lpassedge) mag(rpassedge:length(mag))];
m_stopband = mag(lstopedge:rstopedge);

disp("Max and min passband gain of Kaiser Filter: " + max(m_passband) +", "+ min(m_passband));
disp("Passband variation of Kaiser Filter: " + (max(m_passband)-min(m_passband)))
disp("Peak stopband gain of Kaiser Filter: " + max(m_stopband));

% Equiripple FIR Filter

[n,fo,ao,w] = firpmord(fedges,a,dev,fsamp);
b_pm = firpm(n,fo,ao,w);
disp("Order of Equiripple FIR Filter: " + n);

% Stem plots of Filter coefficients
figure
stem(0:length(b_pm)-1,b_pm,'ro')
title("Equiripple FIR Filter Coefficients")

figure
zplane(b_pm,1)
title('Pole-Zero Plot of Equiripple FIR Filter');
    
% Frequency Response
h = freqz(b_pm,1,f,fsamp);
mag = 20*log10(abs(h));
pha = unwrap(angle(h))*180/pi;
    
figure
subplot(2,1,1)
plot(f/1e6,mag);
xlabel('Frequency (MHz)');
ylabel('Gain (db)');
title('Magnitude Response of Equiripple FIR Filter');
grid on;
axis([0 20 -50 2]);
    
subplot(2,1,2)
plot(f/1e6,pha);
xlabel('Frequency (MHz)');
ylabel('Phase (deg)');
title('Phase Response of Equiripple FIR Filter');
grid on;

% Part c

delratio = delstop/delpass;
wratio = w(1)/w(2); % w(1) gives passband weight, w(2) gives stopband weight
disp("Tolerances have ratio "+ delratio +" and weights for equiripple filter have ratio " + wratio + ": Difference is " + abs(delratio-wratio));
% Seems to check out, difference is within floating-point error

% Part d

m_passband = [mag(1:lpassedge) mag(rpassedge:length(mag))];
m_stopband = mag(lstopedge:rstopedge);

disp("Max and min passband gain of Equiripple FIR Filter: " + max(m_passband) +", "+ min(m_passband));
disp("Passband variation of Equiripple FIR Filter: " + (max(m_passband)-min(m_passband)))
disp("Peak stopband gain of Equiripple FIR Filter: " + max(m_stopband));

% Part e

% The Kaiser window has passband variance much lower than 1.5 dB, as
% achieving the passband variance of 1.5 dB could be done with a lower 
% order filter, but the order was necessary to acheive the desired stopband
% attenuation. 

% Many of the IIR filters also did not acheive 1.5 dB attenuation at the
% "edges", because the left passband edge was at 8.8 MHz and thus would
% expeience much less than 1.5 dB attenuation, while the right passband
% edge was at 12.49 MHz. To mitigate this, I would measure at more points
% and more evenly to get the "edges" closer the the actual frequencies of 9
% and 12.5 MHz.


% Functions

function displayDigitalFilterInfo(ord, freq, h, filtertype)
    lpassedge = find((freq < 9e6),1,"last");
    lstopedge = find((freq > 9.5e6),1,"first");
    rstopedge = find((freq < 12e6),1,"last");
    rpassedge = find((freq > 12.5e6),1,"first");
    disp("Order of "+ filtertype +" Filter: " + ord);
    
    % Frequency Response
    mag = 20*log10(abs(h));
    pha = unwrap(angle(h))*180/pi;
    
    figure
    subplot(2,1,1)
    plot(freq/1e6,mag);
    xlabel('Frequency (MHz)');
    ylabel('Gain (db)');
    title('Magnitude Response of '+ filtertype +' Filter');
    grid on;
    axis([0 20 -50 2]);
    
    subplot(2,1,2)
    plot(freq/1e6,pha);
    xlabel('Frequency (MHz)');
    ylabel('Phase (deg)');
    title('Phase Response of '+ filtertype +' Filter');
    grid on;

    % Peak gain and edge attenuation info
    disp("Peak passband gain of "+ filtertype +" Filter: " + max(mag));
    disp("Attenuation at Stopband Edges of "+ filtertype +" Filter: : " + mag(lstopedge) + " and " + mag(rstopedge));
    disp("Attenuation at Passband Edges of "+ filtertype +" Filter: : " + mag(lpassedge) + " and " + mag(rpassedge));
end