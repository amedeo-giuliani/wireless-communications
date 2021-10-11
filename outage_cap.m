clear;
close all;

epsilon=0.1;
snr_range_db=-10:50;
snr_range=10.^(snr_range_db./10);
fadcap=log2(1+log(1/(1-epsilon)).*snr_range);
awgncap=log2(1+snr_range);

figure;
plot(snr_range_db,fadcap./awgncap);
title('fadcap/awgncap');
hold on;

epsilon=0.01;
fadcap=log2(1+log(1/(1-epsilon)).*snr_range);
plot(snr_range_db,fadcap./awgncap);
hold off;
legend("epsilon=0.1","epsilon=0.01");
grid on;