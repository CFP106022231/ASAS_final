clear; close all;

x = './input/flu.wav';
% [y,fs1] = audioread(x);
% y = y(:,1);
fs = 16000;
% y = resample(y,fs,fs1);
[linAmp, f] = findf0(x, 0.032, 20, 16000);
plot(f)