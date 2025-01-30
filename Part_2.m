clear
clc
%% reading audio
% reading audio
[y,Fs] = audioread('eric.wav');
% fft application
L = length(y); % measure length of function for generation of frequency axis
Y = fftshift(fft(y)); % transform and shift the function
f = (-L/2:L/2-1)*(Fs/L); % generate frequency axis

% plotting the spectrum
plot(f,Y); % plot in frequency domain
title('Frequency Spectrum of eric.wav')
xlabel('Frequency/Hz')
figure;
%% filtering




% filtering
Yfiltered = Y;
Yfiltered(abs(f)>4000) = 0; % simulating ideal lpfilter set at 4000 Hz

% reversing fft and testing
yfiltered = ifft(ifftshift(Yfiltered)); 

% sounding filtered signal
sound(y,Fs); 

pause(8);


%% modulation
% % resampling
yfilresampled = resample(yfiltered,125,12); % resampling to make sampling 
% frequency = carrier frequency times 5

% % modulation
fc = 100000; % set carrier frequency
tcar = 0:1/(fc*5):8.567668; % generate carrier time tupple with overall
% function sampling frequency from 0 to length of audio tupple
tcar(end) = []; % remove final 'extra value'
carrier = cos(2*pi*fc*tcar)'; % generate carrier tupple
 
% % % DSB-SC

sc = yfilresampled.*carrier; % multiply carrier and filtered+resampled function tupple
Ld = length(sc); % calculate legnth of modulated wave
SC = fftshift(fft(sc));
fam = (-Ld/2:Ld/2-1)*((5*fc)/Ld);

plot(fam,SC); % plot DSB-SC spectrum
title('Frequency Spectrum of DSB-SC')
xlabel('Frequency/Hz')
figure

%% Generating LSB

% SC is the DSB-SC spectrum
LSB = SC;
LSB(abs(fam)>fc) = 0; % filter out the USB using ideal lpfilter

plot(fam,LSB);
title('Frequency Spectrum of LSB')
xlabel('Frequency/Hz')
figure

%% Coherent detection

s0 = (ifft(ifftshift(LSB))).*carrier;
S0 = fftshift(fft(s0));
L = length(s0);
f = (-L/2:L/2-1)*(500000/L);
S0(abs(f)>1.5*fc) = 0;
s0 = ifft(ifftshift(S0));

s0R = resample(s0,12,125);
sound(s0R,Fs)
pause(8);
subplot(2,1,1)
plot(tcar,s0)
title('Coherent Detection waveform')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S0)
title('Coherent Detection spectrum')
xlabel('Frequency/Hz')
figure

%% Steps 5 and 6 using Butterworth filter
% step 5 (generating LSB)
[b,a] = butter(4,0.4); % generating Butterworth low pass filter of order 4,
% with Wn = fc/(fc/2) = 100,000/250,000 = 0.4
LSBf = filter(b,a,sc); % generate 'LSB' using Butterworth filter generated
% before

plot(fam,fftshift(fft(LSBf)));
title('Frequency Spectrum of LSB using Butterworth filter')
xlabel('Frequency/Hz')
figure

% step 6 (coherent detection of LSB signal)
sb0 = (ifft(ifftshift(LSBf))).*carrier;
Sb0 = fftshift(fft(sb0));
L = length(sb0);
f = (-L/2:L/2-1)*(500000/L);
Sb0(abs(f)>1.5*fc) = 0;
sb0 = ifft(ifftshift(Sb0));

s0R = resample(s0,12,125);
sound(s0R,Fs)
pause(8);
subplot(2,1,1)
plot(tcar,sb0)
title('Coherent Detection waveform (BW)')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Sb0)
title('Coherent Detection spectrum (BW)')
xlabel('Frequency/Hz')
figure

%% Repeat demodulation with noise

lsbWN0 = awgn(ifft(ifftshift(LSB)),0); % adds white noise 
% with SNR =  0 then  SNR = 10 then SNR = 30
lsbWN10 = awgn(ifft(ifftshift(LSB)),10); 
lsbWN30 = awgn(ifft(ifftshift(LSB)),30);

% with SNR = 0

s0 = lsbWN0.*carrier;
S0 = fftshift(fft(s0));
S0(abs(f)>1.5*fc) = 0;
s0 = ifft(ifftshift(S0));

s0R = resample(s0,12,125);
sound(s0R,Fs)
pause(8);
subplot(2,1,1)
plot(tcar,s0)
title('Coherent Detection waveform SNR = 0')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S0)
title('Coherent Detection spectrum SNR = 0')
xlabel('Frequency/Hz')
figure


% with SNR = 10

s10 = lsbWN10.*carrier;
S10 = fftshift(fft(s10));
S10(abs(f)>1.5*fc) = 0;
s10 = ifft(ifftshift(S10));

s10R = resample(s10,12,125);
sound(s10R,Fs)
pause(8);
subplot(2,1,1)
plot(tcar,s10)
title('Coherent Detection waveform SNR = 10')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S10)
title('Coherent Detection spectrum SNR = 10')
xlabel('Frequency/Hz')
figure

% with SNR = 30

s30 = lsbWN30.*carrier;
S30 = fftshift(fft(s30));
S30(abs(f)>1.5*fc) = 0;
s30 = ifft(ifftshift(S30));

s30R = resample(s30,12,125);
sound(s30R,Fs)
pause(8);
subplot(2,1,1)
plot(tcar,s30)
title('Coherent Detection waveform SNR = 30')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S30)
title('Coherent Detection spectrum SNR = 30')
xlabel('Frequency/Hz')
figure

%% SSB-TC
M = max(yfilresampled); % find amplitude of sound function
tc = (2*M + yfilresampled).*carrier; % generate 'SSB-TC' wave
TC = fftshift(fft(tc));

LSBtc = TC;
LSBtc(abs(fam)>fc) = 0;
lsbtcE = abs(hilbert(ifft(ifftshift(LSBtc))));

plot(tcar,lsbtcE);
title('Waveform of LSB-TC envelope')
xlabel('Time/s')

stcR = resample(lsbtcE,12,125);
sound(stcR,Fs)