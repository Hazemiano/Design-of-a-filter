clc;
clear all;
close all;
% Part 2 %%
%Point 1 : Generate Two Functions for the two linecoding types

function [t, x] = manchester(bits, bitrate)
%encodes BITS array using manchester coding , output is encoded signal
% and time T 

    T = length(bits) / bitrate; % full time of bit sequence
    n = 200;
    N = n * length(bits);
    dt = T / N;
    t = 0:dt:T-dt;
    x = zeros(1, N); % output signal

    for i = 0:length(bits)-1
        if bits(i+1) == 1
            x(i * n+1:(i+0.5) * n) = 1;
            x((i+0.5) * n+1:(i+1) * n) = -1;
        else
            x(i * n+1:(i+0.5) * n) = -1;
            x((i+0.5) * n+1:(i+1) * n) = 1;
        end
    end
end

function [t,x] = unrz(bits, bitrate)
%same as manchester coding function in structure 
% using uni polar non return to zero
T = length(bits)/bitrate; % full time of bit sequence
n = 200;
N = length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t)); % output signal
for i = 0:length(bits)-1
  if bits(i+1) == 1
    x(i * n+1:(i+1) * n) = 1;
  else
    x(i * n+1:(i+1) * n) = 0;
  end
end
end 
% Generate 100k random bits (either 0 or 1)
n_bits = 100000;
random_bits = randi([0, 1], 1, n_bits);

[t1, s1] = manchester(random_bits, 1);
[t2, s2] = unrz(random_bits, 1);

Eb_No_values_db = -10:2:10;
for i = Eb_No_values_db
    Eb_No_linear = 10 ^ (i / 10);
    
    Eb_manchester = 1;
    Eb_unrz = 0.5;
    
    No1 = Eb_manchester / Eb_No_linear;
    No2 = Eb_unrz / Eb_No_linear;
    
    sigma1 = sqrt(No1 / 2 * 200);
    sigma2 = sqrt(No2 / 2 * 200);
    
    noise1 = normrnd(0, sigma1, 1, length(s1));
    noise2 = normrnd(0, sigma2, 1, length(s2));
    
    y1 = s1 + noise1;
    y2 = s2 + noise2;
end
