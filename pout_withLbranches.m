clear;
close all;

L = [1,2,4,8];
R = 1;
snr_range_db = -10:20;
snr_range = 10.^(snr_range_db./10);

pout = zeros(1,length(snr_range_db));

figure;
for l=L
    i = 1;
    for snr=snr_range
        pout(i) = chi2cdf((2^R - 1)/snr,2*l)-chi2cdf(0,2*l);
        i = i+1;
    end
    semilogy(snr_range_db,pout);
    hold on;
end
hold off;
grid on;
legend("L=1","L=2","L=4","L=8");
ylabel("Pout");
xlabel("SNR [dB]");
title("Outage probability vs SNR with diversity L")
ylim([10^-4,1]);