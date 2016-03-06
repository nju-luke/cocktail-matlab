clear;clc;
[wav1,fs]=wavread('mix1.wav');
wav2=wavread('mix2.wav');

X=[wav1 wav2]';

sig=fastica(X,2);

wavwrite(sig(1,:),fs,'sig1.wav');
wavwrite(sig(2,:),fs,'sig2.wav');
