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
pause(9);




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
subplot(2,1,1);

sc = yfilresampled.*carrier; % multiply carrier and filtered+resampled function tupple
Ld = length(sc); % calculate legnth of modulated wave
SC = fftshift(fft(sc)); 
fam = (-Ld/2:Ld/2-1)*((5*fc)/Ld);

plot(fam,SC); % plot DSB-SC spectrum
title('Frequency Spectrum of DSB-SC')
xlabel('Frequency/Hz')

% % % DSB-TC
subplot(2,1,2);

M = max(abs(yfilresampled)); % find amplitude of sound function
tc = (2*M + yfilresampled).*carrier; % generate modulated DSB-TC wave
TC = fftshift(fft(tc));

plot(fam,TC); % plot DSB-TC spectrum
title('Frequency Spectrum of DSB-TC')
xlabel('Frequency/Hz')
figure;

subplot(2,1,1);
% % DSBSC envelope detections
scE = abs(hilbert(sc)); % envelope detection formula
plot(tcar, scE); % plot envelope in time domain
title('DSBSC envelope')
xlabel('Time/s')

subplot(2,1,2);
% % DSBTC envelope detection
tcE = abs(hilbert(tc)); % envelope detection formula
plot(tcar, tcE); % plot envelope in time domain
title('DSBTC envelope')
xlabel('Time/s')
figure

%% resampling and sounding

tcR = resample(tcE,12,125); 
sound(tcR,Fs)
pause(9);
scR = resample(scE,12,125);
sound(scR,Fs)
pause(9);

%% coherent detection for DSB-SC

scWN0 = awgn(sc,0); % adds white noise with SNR = 0 then SNR = 10 then SNR = 30
scWN10 = awgn(sc,10); 
scWN30 = awgn(sc,30);

% % dsbsc with snr = 0
s0 = scWN0.*carrier;
S0 = fftshift(fft(s0));
L = length(s0);
f = (-L/2:L/2-1)*(500000/L);
S0(abs(f)>1.5*fc) = 0;
s0 = ifft(ifftshift(S0));

s0R = resample(s0,12,125);
sound(s0R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,s0)
title('Coherent Detection SNR 0 waveform')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S0)
title('Coherent Detection SNR 0 spectrum')
xlabel('Frequency/Hz')
figure

% % dsbsc with snr = 10
s10 = scWN10.*carrier;
S10 = fftshift(fft(s10));
S10(abs(f)>1.5*fc) = 0;
s10 = ifft(ifftshift(S10));

s10R = resample(s10,12,125);
sound(s10R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,s10)
title('Coherent Detection SNR 10 waveform')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S10)
title('Coherent Detection SNR 10 spectrum')
xlabel('Frequency/Hz')
figure

% % dsbsc with snr = 30
s30 = scWN30.*carrier;
S30 = fftshift(fft(s30));
S30(abs(f)>1.5*fc) = 0;
s30 = ifft(ifftshift(S30));

s30R = resample(s30,12,125);
sound(s30R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,s30)
title('Coherent Detection SNR 30 waveform')
xlabel('Time/s')
subplot(2,1,2)
plot(f,S30)
title('Coherent Detection SNR 30 spectrum')
xlabel('Frequency/Hz')
figure

%% Error 
% % with freq error 100.1 KH
fcf = 100100; % redefining carrier frequency 
carrier = cos(2*pi*fcf*tcar)';
% % dsbsc with snr = 0
se0 = scWN0.*carrier;
Se0 = fftshift(fft(se0));
Se0(abs(f)>1.5*fc) = 0;
se0 = ifft(ifftshift(Se0));

se0R = resample(se0,12,125);
sound(se0R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,se0)
title('Coherent Detection SNR 0 waveform and freq error')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Se0)
title('Coherent Detection SNR 0 spectrum and freq error')
xlabel('Frequency/Hz')
figure

% error calculation for snr = 0
errorArr0 = se0/s0;
plot(tcar,errorArr0);
title('Freq error for SNR = 0')
xlabel('Time/s')
figure

% % dsbsc with snr = 10
se10 = scWN10.*carrier;
Se10 = fftshift(fft(se10));
Se10(abs(f)>1.5*fc) = 0;
se10 = ifft(ifftshift(Se10));

se10R = resample(se10,12,125);
sound(se10R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,se10)
title('Coherent Detection SNR 10 waveform and freq error')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Se10)
title('Coherent Detection SNR 10 spectrum and freq error')
xlabel('Frequency/Hz')
figure

% error calculation for snr = 10
errorArr10 = se10/s10;
plot(tcar,errorArr10);
title('Freq error for SNR = 10')
xlabel('Time/s')
figure

% % dsbsc with snr = 30
se30 = scWN30.*carrier;
Se30 = fftshift(fft(se30));
Se30(abs(f)>1.5*fc) = 0;
se30 = ifft(ifftshift(Se30));

se30R = resample(se30,12,125);
sound(se30R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,se30)
title('Coherent Detection SNR 30 waveform with freq error')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Se30)
title('Coherent Detection SNR 30 spectrum with freq error')
xlabel('Frequency/Hz')
figure

% error calculation for snr = 30
errorArr30 = se30/s30;
plot(tcar,errorArr30);
title('Freq error for SNR = 30')
xlabel('Time/s')
figure

%% Error 
% % with phase error 20
carrier = cos(2*pi*fc*tcar + 20)'; % generating carrier again with phase error
% % dsbsc with snr = 0
sp0 = scWN0.*carrier;
Sp0 = fftshift(fft(sp0));
L = length(sp0);
f = (-L/2:L/2-1)*(500000/L);
Sp0(abs(f)>1.5*fc) = 0;
sp0 = ifft(ifftshift(Sp0));

sp0R = resample(sp0,12,125);
sound(sp0R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,sp0)
title('Coherent Detection SNR 0 waveform with phase error')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Sp0)
title('Coherent Detection SNR 0 spectrum with phase error')
xlabel('Frequency/Hz')
figure


% % dsbsc with snr = 10
sp10 = scWN10.*carrier;
Sp10 = fftshift(fft(sp10));
L = length(sp10);
f = (-L/2:L/2-1)*(500000/L);
Sp10(abs(f)>1.5*fc) = 0;
sp10 = ifft(ifftshift(Sp10));

sp10R = resample(sp10,12,125);
sound(sp10R,Fs)
pause(9);

subplot(2,1,1)
plot(tcar,sp10)
title('Coherent Detection SNR 10 waveform with phase error')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Sp10)
title('Coherent Detection SNR 10 spectrum with phase error')
xlabel('Frequency/Hz')
figure

% % dsbsc with snr = 30
sp30 = scWN30.*carrier;
Sp30 = fftshift(fft(sp30));
L = length(sp30);
f = (-L/2:L/2-1)*(500000/L);
Sp30(abs(f)>1.5*fc) = 0;
sp30 = ifft(ifftshift(Sp30));

sp30R = resample(sp30,12,125);
sound(sp30R,Fs)

subplot(2,1,1)
plot(tcar,sp30)
title('Coherent Detection SNR 30 waveform with phase error')
xlabel('Time/s')
subplot(2,1,2)
plot(f,Sp30)
title('Coherent Detection SNR 30 spectrum with phase error')
xlabel('Frequency/Hz')