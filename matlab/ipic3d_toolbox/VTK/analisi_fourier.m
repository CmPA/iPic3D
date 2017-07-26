function [f, out, NFFT]=analisi_fourier(y, Dt)

Fs = 1/Dt;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = max(size(y));             % Length of signal
t = (0:L-1)*T;                % Time vector

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);


out = Y(1:round(NFFT/2+1));

end


