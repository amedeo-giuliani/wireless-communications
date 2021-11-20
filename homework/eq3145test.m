clear;
close all;

%to test relationship (3.145) of the book

W = 2*10^6; %bandwidth
fc = 900*10^6; %carrier frequency
psksize = 4;
h = [0.7; 0.25; 0; 0; 0.4; 0.5; 0; 0.1; 0.1];
h = h./norm(h); %normalize channel impulse response such that sum(h)=1
L = length(h);
CP = L-1;
N = 2^ceil(log2(9*CP)); %choose N such as CP/(CP+N) <= 10%
n = 0:N-1;
f = fc + n*W/N;
figure;
stem(h);
title('Channel impulse response long L=5');
figure;
H = fft(h,N);
stem(f,abs(H));
ylim([0 inf]);
xlabel('[Hz]');
title('Modulus of channel frequency response vs frequency');
pause;

SNR_range_dB = 0:35;
%split the total power equally among all the subcarriers as TX does not have CSIT(not optimal, of course)
P_opt = ones(N,1); %assume unitary total tx power -> 1/N power to each subcarrier

bervec = zeros(1,length(SNR_range_dB));
spctreff = zeros(1,length(SNR_range_dB));

for SNR_dB = SNR_range_dB
    d0 = randi([0 1],N,1);
%     d = pskmod(d0,2);
    D = sqrt(N) * ifft(d0,N); %go back in the time domain; sqrt(N) is to preserve transmit power
    x = [D(end-CP+1:end); D]; %add the cyclic prefix
    y = filter(h,1,x);
    y = reshape(y,N+CP,[]);
    R = y(CP+1:end);
    r0 = fft(R,N)/sqrt(N);
    r_zf = r0./H;
%     r = pskdemod(r_zf,2);
    r = real(r_zf);
    r = cast(r,'uint8');
%     ber = biterr(d0,r_zf);
    ber = biterr(d0,r);
end

figure;
stem(real(r0),'or');
hold on
q = H.*d0;
stem(real(q),'xk');
hold off
legend('rx symbols','H.*d0 symbols');
title('real part');

figure;
stem(imag(r0),'or');
hold on
stem(imag(q),'xk');
hold off
legend('rx symbols','H.*d0 symbols');
title('imaginary part');