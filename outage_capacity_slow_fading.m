clear;
close all;

snr_range = -10:40; %dB
epsilon = 0.1;
fadcap = zeros(1,length(snr_range));
Q = 10000; %average pout over Q realizations of h
R = zeros(1,length(snr_range));
%to speed up convergence
R(21:29)=1;
R(30:38)=3;
R(39:41)=6;
R(42:44)=7;
R(45:48)=8;
R(48:end)=9;
snr = 10.^(snr_range./10); %linear

while true
    pout = zeros(1,length(snr_range));
    for i=1:Q
        h = 1/sqrt(2)*(randn()+1i*randn()); %rayleigh flat slow fading ch
        idx = find(log2(1+abs(h)^2*snr) < R);
        pout(idx) = pout(idx)+1; %outage if info rate is larger than actual capacity
    end
    pout=pout./Q;
    idx = find(pout < epsilon-0.02); %to make the solution closer to epsilon
    if size(idx,2)==0 %if all the pouts are below epsilon, we're done
        break;
    end
    R(idx)=R(idx)+0.001;
end

disp(pout); %outage probability versus SNR which we made it to be epsilon
disp(R); %epsilon outage capacity versus SNR
awgncap=log2(1+snr);
plot(snr_range,R./awgncap); % epsilon outage capacity over awgn capacity

% for q=1:Q
%     i=1;
%     h = 1/sqrt(2)*(randn()+1i*randn()); %flat slow fading ch
%     for snrdb=snr_range
%         snr = 10^(snrdb/10);
%         fadcap(i) = fadcap(i)+log2(1+abs(h)^2*snr);
%         awgncap(i) = awgncap(i)+log2(1+snr);
%         i=i+1;
%     end
% end
% 
% fadcap=fadcap./Q;
% awgncap=awgncap./Q;
% 
% figure;
% plot(snr_range,fadcap./awgncap);
% title('fadcap/awgncap');