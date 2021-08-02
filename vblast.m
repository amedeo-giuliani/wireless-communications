clear;
close all;

M = 4; %modulation order
Ntx = 8; %number of tx antennas
Nrx = 8; %number of rx antennas
P_tx = 1; %max tx power per antenna
K = 1e2; %number of independent experiments
SNR_rangedB = -10:20;
cap = zeros(1,length(SNR_rangedB));
servec = zeros(1,length(SNR_rangedB));

for k = 1:K
    idx = 1;
    i = 1;
    for SNR_dB = SNR_rangedB
        H = 1/sqrt(2)*(randn(Nrx,Ntx) + 1i*randn(Nrx, Ntx)); %random iid rayleigh MIMO
        [~,L,~] = svd(H);
        No = 10^(-SNR_dB/10);
        P_opt = Ntx;
        if Ntx > 1
            %P_opt = waterfilling(diag(L).^2,No,Ntx*P_tx);
            P_opt = ones(1,Ntx)/Ntx;
        end
        x0 = randi([0 3],Ntx,1); %generate quaternary symbol for each tx antenna
        x = qammod(x0,M); %M-QAM modulation
        y = awgn(H*x,SNR_dB); %fading and additive ZMCSCG noise
        
        %maximal ratio combiner receiver when only one tx antenna
        if Ntx == 1
            mrc = H'./sqrt(sum(abs(H').^2));
            r = qamdemod(mrc*y,M);
        end
        
        %zero-forcing sic receiver as in vblast literature
        if Ntx > 1
            H1 = H;
            r = zeros(Nrx, 1);
            orderVec = 1:Nrx;
            j = Nrx+1;
            % Start ZF nulling loop
            for n = 1:Nrx
                % Shrink H to remove the effect of the last decoded symbol
                H1 = H1(:, [1:j-1,j+1:end]);
                % Shrink order vector correspondingly
                orderVec = orderVec(1, [1:j-1,j+1:end]);
                % Select the next symbol to be decoded
                G = (H1'*H1) \ eye(Nrx-n+1); % Same as inv(H'*H), but faster
                [~, j] = min(diag(G));
                symNum = orderVec(j);
                
                % Hard decode the selected symbol
                dec = qamdemod(G(j,:) * H1' * y,M);
                r(symNum) = dec;
                
                % Subtract from y the effect of the last decoded symbol
                if n < Nrx
                    y = y - H1(:, j) * qammod(dec,M);
                end
            end
        end
        [~,ser] = symerr(x0,r);
        servec(idx) = servec(idx) + ser;
        
        %MIMO capacity formula
        if Ntx > 1
            cap(idx) = cap(idx) + sum(log2(1 + P_opt.'.*diag(L).^2/No));
        end
        
        %SIMO capacity formula
        if Ntx == 1
            cap(idx) = cap(idx) + log2(1 + P_opt * abs(sum(H)).^2 / No);
        end
        idx = idx + 1;
    end
end

servec = servec./K;
cap = cap./K;

figure;
semilogy(SNR_rangedB,servec,'-xr');
xlabel('SNR [dB]');
ylabel('SER');
title('SER vs SNR');
grid on;

figure;
plot(SNR_rangedB, cap, '+g-');
xlabel('SNR [dB]');
%ylim([0, 1]);
ylabel('Capacity [bit/s/Hz]');
title('Capacity vs SNR');
grid on;

function Popt = waterfilling(H2, N0, barP)


N = length(H2);

[H2sort, idx] = sort(H2,'descend'); % order channel gains in descending order


for n=N:-1:1
    
    lambda = n/(barP + sum(N0./H2sort(1:n)));
    P = (1/lambda - N0./H2sort(1:n)); % compute powers
    
    if (P(end)>=0) %found solution
        Pbest=P;
        pos=n;
        break
    end
    % otherwise consider less subcarriers
end

Popt = zeros(1,N);
Popt(idx(1:n))= Pbest;
end