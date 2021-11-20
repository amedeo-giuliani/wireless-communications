clear all;
close all;
warning('off','all');

W = 2*10^6; %bandwidth
fc = 900*10^6; %carrier frequency
psksize = 4;
h = [0.7; 0.25; 0; 0; 0.4; 0.5; 0; 0.1; 0.1];
h = h./norm(h); %normalize channel impulse response such that sum(h)=1
L = length(h);
CP = L-1;
N = 2^ceil(log2(9*CP)); %choose N such as CP/(CP+N) <= 10%
n = -N/2:N/2-1;
f = fc + n*W/N;
figure;
stem(h);
title('Channel impulse response');
figure;
H = fft(h,N);
stem(f./1e6,abs(H));
ylim([0 inf]);
xlabel('[MHz]');
title('Magnitude of channel frequency response vs frequency');
m = 1000; %number of ofdm blocks to transmit
SNR_range_dB = 0:30;
%split the total power equally among all the subcarriers as TX does not have CSIT(not optimal, of course)
P_opt = ones(N,1); %assume total tx power equal to N -> unitary power to each subcarrier

bervec = zeros(1,length(SNR_range_dB));
spctreff = zeros(1,length(SNR_range_dB));
ub = zeros(1,length(SNR_range_dB));
lb = zeros(1,length(SNR_range_dB));
rxbits=zeros(1,length(SNR_range_dB));
% if we apply waterfilling algo, then we can use a higher order modulation
% on the subcarriers with high SNR and lower order modulation on the ones
% with low SNR, i.e. we transmit more on the good subcarriers

qpskmod = comm.QPSKModulator('BitInput',true,'SymbolMapping','Gray'); %use QPSK to modulate the bits with Gray mapping to minimize the distance between symbols -> maximize BER
qpskdemod = comm.QPSKDemodulator('BitOutput',true,'SymbolMapping','Gray'); %demodulate QPSK symbols

%modulation orders to modulate bits after waterfilling
q1 = 64;
q2 = 16;
q3 = 16;
q4 = 4;

%flag to use waterfilling or not
wf = 1;

K = 10;
for k=1:K
    disp(k);
    i = 1;
    for SNR_dB = SNR_range_dB
        No = 10^(-SNR_dB/10); %noise power at current SNR
        if wf==1
            P_opt = waterfilling(abs(H).^2,No,N)';
            maxP=max(P_opt);
            firstmod = find(P_opt>=maxP-maxP/4);
            secondmod = find(P_opt<maxP-maxP/4 & P_opt>=maxP-maxP/2.1);
            %thirdmod = find(P_opt<maxP-maxP/3 & P_opt>=maxP-maxP/2);
            fourthmod = find(P_opt<maxP-maxP/2.1 & P_opt>=0);
            %use a constellation scaled with the average power allocated to
            %the subcarriers that belong to firstmod,secondmod,thirdmod and
            %fourthmod
            mod1 = comm.RectangularQAMModulator(q1,'BitInput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(firstmod)));
            mod2 = comm.RectangularQAMModulator(q2,'BitInput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(secondmod)));
%             mod3 = comm.RectangularQAMModulator(q3,'BitInput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(thirdmod)));
            mod4 = comm.RectangularQAMModulator(q4,'BitInput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(fourthmod)));
        end
        for j=1:m
            if wf==0
                bits = randi([0 1],log2(psksize)*N,1); %generate bits
                d = qpskmod(bits);
            end
            if wf==1
                bits = randi([0 1],log2(q1),N); % generate at most one matrix of log2(highest mod order) x N
                a1 = bits(:,firstmod);
                a2 = bits(:,secondmod);
                %a3 = bits(:,thirdmod);
                a4 = bits(:,fourthmod);
                d1=[];
                d2=[];
                %d3=[];
                d4=[];
                for t=1:size(a1,2)
                    d1(end+1)= mod1(a1(1:log2(q1),t));
                end
                for t=1:size(a2,2)
                    d2(end+1)= mod2(a2(1:log2(q2),t));
                end
%                 for t=1:size(a3,2)
%                     d3(end+1)= mod3(a3(1:log2(q3),t));
%                 end
                for t=1:size(a4,2)
                    d4(end+1)= mod4(a4(1:log2(q4),t));
                end
                d = zeros(1,N); % prepare the vector to transmit
                % fill it up with all the modulated symbols mapped back to the
                % subcarrier they belong to
                d(firstmod) = d1;
                d(secondmod) = d2;
                %d(thirdmod) = d3;
                d(fourthmod) = d4;
                d = reshape(d,N,1); % reshape the row vector into a column vector
            end
            D = sqrt(N) * ifft(d,N); %go back in the time domain; *sqrt(N) is to preserve the power
            x = [D(end-CP+1:end,:); D]; %add the cyclic prefix
            
            y = filter(h,1,x); %ofdm block goes through the fading channel
            y = awgn(y,SNR_dB); %then it goes through the awgn channel
            R = y(CP+1:end,:); %remove the cyclic prefix
            r0 = fft(R,N)/sqrt(N); %go back in frequency domain; /sqrt(N) is to preserve the power
            r_zf = r0./H; %zero-forcing receiver
            if wf==0
                r = qpskdemod(r_zf);
                err = comm.ErrorRate();
                e = err(bits,r);
                ber = e(1);
                numerr = e(2);

                %compute CIs at level 95%(of course they'll be larger and larger with
                %growing SNR as the number of errors becomes smaller)
                [~,ci] = berconfint(numerr,length(bits),0.95);
                lb(i) = lb(i) + ci(1);
                ub(i) = ub(i) + ci(2);
                rxbits(i) = rxbits(i) + length(bits);
            end
            
            if wf==1
                demod1 = comm.RectangularQAMDemodulator(q1,'BitOutput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(firstmod)));
                demod2 = comm.RectangularQAMDemodulator(q2,'BitOutput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(secondmod)));
                %demod3 = comm.RectangularQAMDemodulator(q3,'BitOutput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(thirdmod)));
                demod4 = comm.RectangularQAMDemodulator(q4,'BitOutput',true,'NormalizationMethod','Average Power','AveragePower',mean(P_opt(fourthmod)));
                rx1 = r_zf(firstmod);
                rx2 = r_zf(secondmod);
                %rx3 = r_zf(thirdmod);
                rx4 = r_zf(fourthmod);
                r1=[];
                r2=[];
                %r3=[];
                r4=[];
                for t=1:length(rx1)
                    r1(:,end+1)=demod1(rx1(t));
                end
                for t=1:length(rx2)
                    r2(:,end+1)=demod2(rx2(t));
                end
%                 for t=1:length(rx3)
%                     r3(:,end+1)=demod3(rx3(t));
%                 end
                for t=1:length(rx4)
                    r4(:,end+1)=demod4(rx4(t));
                end
                [~,s1] = size(a1);
                [~,s2] = size(a2);
%                 [~,s3] = size(a3);
                [~,s4] = size(a4);
                e1=0;
                e2=0;
%                 e3=0;
                e4=0;
                if s1>0
                    [e1,ber1] = biterr(a1,r1);
                end
                if s2>0
                    [e2,ber2] = biterr(a2(1:log2(q2),:),r2);
                end
%                 if s3>0
%                     [e3,ber3] = biterr(a3(1:log2(q3),:),r3);
%                 end
                if s4>0
                    [e4,ber4] = biterr(a4(1:log2(q4),:),r4);
                end
                %ber = (e1+e2+e3+e4)/(numel(a1)+numel(a2(1:log2(q2),:))+numel(a3(1:log2(q3),:))+numel(a4(1:log2(q4),:)));
%                 ber = (e1+e2+e3+e4)/(numel(a1)+numel(r2)+numel(r3)+numel(r4));
                ber = (e1+e2+e4)/(numel(a1)+numel(r2)+numel(r4));
                rxbits(i) = rxbits(i) + numel(r1)+numel(r2)+numel(r4);
            end
            
            bervec(i) = bervec(i) + ber;
            spctreff(i) = spctreff(i) + sum(log2(1+abs(H).^2.*P_opt./No))/(N+CP);
        end
        i = i + 1;
    end
end

bervec = bervec./(m*K);
spctreff = spctreff./(m*K);
capacity = spctreff.*W;
throughput = rxbits./(m*K);
ub = ub./(m*K);
lb = lb./(m*K);

if wf==0
    figure;
    semilogy(SNR_range_dB,bervec,'-xr');
    hold on
    errorbar(SNR_range_dB,bervec,bervec-lb,ub-bervec,'LineStyle','none','Color','k');
    hold off
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER vs SNR with CIs at level 95%');
    grid on;
end
if wf==1
    figure;
    semilogy(SNR_range_dB,bervec,'-xr');
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER vs SNR');
    grid on;
end

figure;
plot(SNR_range_dB, capacity, '+g-');
xlabel('SNR (dB)');
%ylim([0, 1]);
ylabel('Capacity [bit/s]');
title('Capacity vs SNR');
grid on;

figure;
plot(SNR_range_dB, throughput, '+g-');
xlabel('SNR (dB)');
%ylim([0, 1]);
ylabel('Throughput [bit/s]');
title('Throughput vs SNR');
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