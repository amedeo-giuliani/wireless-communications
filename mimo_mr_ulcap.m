clear;
close all;

K = 10; %number of users
beta = 1; %large scale fading factor
tau = K; %one pilot per user
eta = 1; %ul power control for each user
M = 1:500; %number of antennas
ro_ul_dB = [-10,-5,0,5]; %uplink power in dB for each user
ro_ul = 10.^(-ro_ul_dB/10);
gamma = tau.*ro_ul.*beta^2 ./ (1 + tau.*ro_ul.*beta);

C_mr = zeros(length(M),length(ro_ul)); %max rate for each user

i = 1;
for m=M
    j = 1;
    for r=ro_ul
        C_mr(i,j) = log2(1 + m*r*eta*gamma(find(ro_ul==r)) / (1 + sum(repelem(r*eta*beta,K))));
        j = j+1;
    end
    i = i+1;
end

figure;
for C=C_mr
    plot(M,C);
    hold on;
end
hold off;
grid on;
legend('ro = 5 dB','ro = 0 dB','ro = -5 dB','ro = -10 dB','Location','northwest');
title('Uplink MR capacity lower bound with varying SNR');
ylabel('Capacity [bit/s/Hz]');
xlabel('Number of antennas');

K = 1:20;
tau = K;
M1 = 100;
ro_ul_dB = 5;
ro_ul = 10^(-ro_ul_dB/10);
gamma = tau.*ro_ul.*beta^2 ./ (1 + tau.*ro_ul.*beta);

C_mr = zeros(length(K),2);

i = 1;
for k=K
    M2 = 10.*k;
    C_mr(i,1) = log2(1 + M1*ro_ul*eta*gamma(find(tau==k)) / (1 + sum(repelem(ro_ul*eta*beta,k))));
    C_mr(i,2) = log2(1 + M2*ro_ul*eta*gamma(find(tau==k)) / (1 + sum(repelem(ro_ul*eta*beta,k))));
    i = i+1;
end

figure;
for C=C_mr
    plot(K,C);
    hold on;
end
hold off;
grid on;
legend('M = 100','M = 10K','Location','southeast');
title('Uplink MR capacity lower bound with scaling number of users');
ylabel('Capacity [bit/s/Hz]');
xlabel('Number of users');

C_mr_sum = zeros(length(K),2);

i = 1;
for k=K
    C_mr_sum(i,:) = k.*C_mr(i,:);
    i = i+1;
end

figure;
for C=C_mr_sum
    plot(K,C);
    hold on;
end
hold off;
grid on;
legend('M = 100','M = 10K','Location','northwest');
title('Uplink MR sum capacity lower bound with scaling number of users');
ylabel('Capacity [bit/s/Hz]');
xlabel('Number of users');