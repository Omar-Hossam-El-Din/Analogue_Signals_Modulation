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
% % sound(y,Fs);
% % pause(8);
%% modulation
% by Narrow-band FM
% % resampling
yfilresampled = resample(yfiltered,125,12); % resampling to make sampling 
% frequency = carrier frequency times 5


% % modulation
fc = 100000; % set carrier frequency
tcar = 0:1/(fc*5):8.567668; % generate carrier time tupple with overall
% function sampling frequency from 0 to length of audio tupple
tcar(end) = []; % remove final 'extra value'
Ld = length(yfilresampled); % calculate length of resampled wave
ffm = (-Ld/2:Ld/2-1)*((5*fc)/Ld); 

kf = 0.1; % set frequency variation constant to low enough value
x = (cumsum(yfilresampled))'; % calculate cumulative summation of audio 
% function from -inf to 't'
s = cos(2*pi*fc*tcar + 2*pi*kf*x); % generate NBFM signal
S = fftshift(fft(s)); % fourier transform and shift


plot(ffm,S); % plot FM modulated signal in frequency domain
title('Frequency Spectrum of FM modulated signal')
xlabel('Frequency/Hz')
figure

%% demodulation using differentiator and envelope detector

sdiff = diff(s);
senv = abs(hilbert(sdiff)); % envelope detector function
tcar(end) = [];
plot(tcar,senv)
title('Waveform of Envelope')
xlabel('Time/s')

