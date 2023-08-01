%this code generates a combination of two sinusoidal inputs and
%write them into a .dat file in binary format

clc;
clear all;
close all;

fs = 2*10^3;  %sampling frequency = 1kHz
f1 = 0.2*10^3;  %input1 frequency = 200Hz
f2 = 0.6*10^3;  %input2 frequency = 600Hz

size = 16; %16 bit input
sf = 2^15; %scaling factor

n = 1:fs;
x_n1 = 0.5*sin((2*pi*f1*n)/fs);
x_n2 = 0.5*sin((2*pi*f2*n)/fs);
x_n = x_n1+x_n2;


FFT_xn = fftshift(fft(x_n)); %FFT for the input signals
freq = -fs/2:fs/2000:(fs/2-fs/2000);

figure();
plot(freq,abs(FFT_xn)); %plotting the FFT
title('FFT of Input');
xlabel('Frequency(Hz)');
ylabel('Magnitude');

%conversion of the input samples to binary 16-bit format
x_n_int = round(x_n.*sf);
x_n_binary = dec2bin(mod((x_n_int),2^size),size);
x_n_b = cellstr(x_n_binary);

%writing the binary 16bit output in a file
fid = fopen('input_comb.dat','wt');
fprintf(fid, '%8s \n', x_n_b{:});
disp('text file for signal finished');