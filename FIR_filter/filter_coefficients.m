%This code will read the filter coefficients and convert them to signed
%binary format.

clc;
clear all;
close all;

size = 16; %size of the word
sf = 2^15; %scaling factor

load('filter_coeff.mat');

coeff_integers = round(filter_coeff.*sf);
coeff_binary = dec2bin(mod((coeff_integers),2^size),size);
coeff_b = cellstr(coeff_binary);

fid = fopen('coeff.dat','wt');
fprintf(fid, '%8s \n', coeff_b{:});
disp('text file for coefficients finished');