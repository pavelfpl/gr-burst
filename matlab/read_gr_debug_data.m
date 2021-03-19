clear
clc

close all
% Read gr debug data ...
% ----------------------
% fname = '/tmp/gr_inputData.txt';
fname = '/tmp/gr_inputDataDecimated.txt';

inputData = dlmread(fname);
inputData = inputData(1:2:end) + 1.0j*inputData(2:2:end);
inputData = inputData.';

% load good_burst
% inputData = x;

Fs = 97.65625e3;
debugFilename = 0;

[pll_v1_out, pll_v2_out, pll_v3_out,pll_v4_out, preCrossCorr] = ...
        qpskSyncBurst(inputData, Fs, .002, debugFilename);
 
% return    
    
figure;
subplot(2,3,1);
plot(preCrossCorr)
grid on
title('Burst Preamble Cross Correlation')

subplot(2,3,2);
scatter(real(pll_v1_out),imag(pll_v1_out))
grid on
title('Constellation w/ v1 PLL')

subplot(2,3,3);
scatter(real(pll_v2_out),imag(pll_v2_out))
grid on
title('Constellation w/ v2 PLL')

subplot(2,3,4)
scatter(real(pll_v3_out),imag(pll_v3_out))
grid on
title('Constellation w/ v3 PLL')

subplot(2,3,5)
scatter(real(pll_v4_out),imag(pll_v4_out))
grid on
title('Constellation w/ v4 (second order) PLL')