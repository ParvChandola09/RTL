clc;
clear all;
close all;

% Define the length of the array
arrayLength = 32;
size = 16; %size of the word
sf = 2^8; %scaling factor

% Generate a random array with values between 0 and 10
input_dec = randi([0, 5], 1, arrayLength);

input_scaled = input_dec.*sf;
bin_input = dec2bin(mod((input_scaled),2^size),size);
input_binary = cellstr(bin_input);

filename = 'input_decimal.dat';

% Open the file for writing
fileID = fopen(filename, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Could not open the file for writing.');
end

% Loop through the array and write each value to the file
for i = 1:numel(input_dec)
    fprintf(fileID, '%d\n', input_dec(i));
end

% Close the file
fclose(fileID);

%writing the binary input
% Specify the filename for the .dat file
filename = 'input_data.dat';

% Open the file for writing
fileID = fopen(filename, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Could not open the file for writing.');
end

% Loop through the cell array and write each binary string to the file
for i = 1:numel(input_binary)
    fprintf(fileID, '%s\n', input_binary{i});
end

% Close the file
fclose(fileID);

