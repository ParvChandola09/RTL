clc;
clear all;
close all;

sf = 2^-8;

%reading the decimal input
data_decimal = importdata('input_decimal.dat');

%reading the RTL output
data_RTL = importdata('FFT32.txt') ; 
D = data_RTL.*sf;
RTL_output = complex(D(:,1),D(:,2));

%MATLAB output
matlab_output = fft(data_decimal);

figure();
subplot(1,2,1);
stem(abs(matlab_output));
title('MATLAB output');
grid on;
subplot(1,2,2);
stem(abs(RTL_output));
title('RTL output');
grid on;
