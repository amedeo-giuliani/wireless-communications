clear;
close all;

K = [1,2,4,16]; % number of users in the cell
Ntx = 1; % number of antennas at BS
Nrx = 1; % number of antennas at UE
P = 1; % power assigned to each user
snr_range_db = -20:20;

C_sumCSIR = zeros(length(snr_range_db),length(K));
C_sumCSIT = zeros(length(snr_range_db),length(K));

j = 1;
for k=K
    i = 1;
    for snr_db=snr_range_db
        No = 10^(-snr_db/10);
        for t = 1:1e4
            H = 1/sqrt(2)*(randn(Nrx*k,Ntx) + 1i*randn(Nrx*k, Ntx)); % fast fading rayleigh channel between BS and users
            C_sumCSIR(i,j) = C_sumCSIR(i,j) + log2(1 + P*sum(abs(H).^2)/No);
            
            [~,idx] = max(H);
            C_sumCSIT(i,j) = C_sumCSIT(i,j) + log2(1 + k*P*abs(H(idx))^2/No);
        end
        i = i + 1;
    end
    j = j + 1;
end

C_sumCSIR = C_sumCSIR ./ t;
C_sumCSIT = C_sumCSIT ./ t;

% C_sumCSIT = zeros(length(snr_range_db),length(K));
% 
% j = 1;
% for k=K
%     i = 1;
%     for snr_db=snr_range_db
%         No = 10^(-snr_db/10);
%         for t = 1:1e5
%             H = 1/sqrt(2)*(randn(Nrx*k,Ntx) + 1i*randn(Nrx*k, Ntx)); % fast fading rayleigh channel between BS and users
%             [~,idx] = max(H);
%             C_sumCSIT(i,j) = C_sumCSIT(i,j) + log2(1 + k*P*abs(H(idx))^2/No);
%         end
%         i = i + 1;
%     end
%     j = j + 1;
% end
% 
% C_sumCSIT = C_sumCSIT ./ t;

C_awgn = zeros(length(snr_range_db),length(K));
j = 1;
for k=K
    i = 1;
    for snr_db=snr_range_db
        No = 10^(-snr_db/10);
        C_awgn(i,j) = log2(1 + k*P/No);
        i = i + 1;
    end
    j = j + 1;
end

figure;
plot(snr_range_db,C_sumCSIR(:,1),'--xb');
hold on;
plot(snr_range_db,C_awgn(:,4),'--og');
hold on;
plot(snr_range_db,C_sumCSIT(:,1),'--+r');
hold on;
legend('CSIR','AWGN','CSIT','AutoUpdate','off');
plot(snr_range_db,C_sumCSIR(:,2),'--xb');
hold on;
plot(snr_range_db,C_sumCSIR(:,3),'--xb');
hold on;
plot(snr_range_db,C_sumCSIR(:,4),'--xb');
hold on;
% plot(snr_range_db,C_awgn(:,2),'--g');
% hold on;
% plot(snr_range_db,C_awgn(:,3),'--g');
% hold on;
% plot(snr_range_db,C_awgn(:,4),'--g');
% hold on;
plot(snr_range_db,C_sumCSIT(:,2),'--+r');
hold on;
plot(snr_range_db,C_sumCSIT(:,3),'--+r');
hold on;
plot(snr_range_db,C_sumCSIT(:,4),'--+r');
grid('on');