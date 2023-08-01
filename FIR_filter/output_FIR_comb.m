%This code will plot the FFT of the output

clc;
clear all;
close all;

sf = 2^-15;
fs = 2000;

data = importdata('filter_output_comb1.txt') ; 
output_original = data.*sf;

FFT_xn = fftshift(fft(output_original));
freq = -fs/2:fs/2000:(fs/2-fs/2000);

figure();
plot(freq,abs(FFT_xn)); %plotting the FFT
title('FFT of Output');
xlabel('frequency(Hz)');
ylabel('Magnitude');