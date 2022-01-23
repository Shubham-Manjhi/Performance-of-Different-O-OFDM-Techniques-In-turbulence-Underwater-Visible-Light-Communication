%%%%%%% OFDM-ACO
clear all
clc
N=64; %number of subcarrier
%  generating random data symbols
data_bits= randi([0,3],64000,1);

%% Parameters
%using FEC or not
FEC_Coding = 0;

%SNR 
N0 = 10^-21;  % spectral density level of gaussian noise at the receiver
B = 20*10^12;  % signal bandwidth B = 20MHz

SNR = 100:5:140; %dB, Desired electrical SNR per bit before clipping at the transmitter (excluding the dc_bias)


%ACO-OFDM params
Nsym = 10000;  %number of OFDM symbols
subcar = 32;
FFT_size = subcar*4;
cp_size = 16;

nBitPerSymbol = 4;
yita = nBitPerSymbol*subcar/(4*subcar+cp_size); % Spectral efficiency

%LED params
led.min = 0.1; % A
led.max = 1 ; % A
led.dc_bias = led.min; %set the dc bias according to the signal variance

%LED filter (Tx filter)
% led_filter = [];
% %VLC channel filter
% vlc_filter = [];
% channel = conv(led_filter,vlc_filter);
channel = 4*10^-6; %Dirac channel

%Rx filter
pd = 1;


% QAM symbol mapping
 symbols= qammod(data_bits,4);
%
x_mod=reshape(symbols,N/4,length(symbols)/(N/4));
snr=0:30;
for t=1:length(snr)
    
    
no_of_error=[];
ratio=[];
%TO1=toc(Tin1);   
    
    for i=1:(length(symbols)/(N/4))
        
        
        
        x1=x_mod(:,i);
        x_herm=flipud(conj(x1)); 
        x_temp= [x1; x_herm];
        x_zero=[0;upsample(x_temp,2)];
         x_final= x_zero(1:end-1); 
         
        
         
        x_ifft=ifft(x_final,N); 

        %CP ADDITION
          x_ifft= [x_ifft(1 : N/4); x_ifft];
          
        %clipping
        %%x_ifft(x_ifft<0)=0; 
        %%x_ifft(x_ifft>1)=1; 
        %%
        
        %% LED Clipping  
        led_clipping=@aco_ofdm_led_filter;
        [txData,ChTrainModel] = vlc_led_filter(x_ifft, led_clipping, led);
        x11(:,i)= x_ifft;
        
        
        
        %% VLC Channel

          alpha = 10;% beam half angle in degree
          psi =pi*(alpha*pi/180)^2; % alpha is converted first in radian and then calculate beam solid angle i.e. pi*(alpha)^2
          mode=45; %mode number
          Ar = 9.8*10^(-6);% receiver PD area= 9.8mm^2
          R = 1;% separation= 5000 mm
          h= ((mode+1)/2*pi)*(Ar/R^2);

          r= 1;
          packet = r*(h*x_ifft);
          rx = awgn(packet,snr(t),'measured');

          
          
          
        %% Data through receiver
        Rx=rx/h;
        Rx= Rx(N/4+1:end);
        y_fft=fft(Rx,N);
        % y_fft=y_fft(1:length(x)/k+1);
         y_fft=y_fft(2:end);
        % y_fft= [y_fft(2:end);0];
        y_final= downsample(y_fft,2);
        Len = length(y_final);
        y_final= y_final(1:Len/2);

        %% Demodulation
        y_demod(:,i) = qamdemod(y_final,4);
end

    
    
rec_symbols = reshape(y_demod,length(symbols),1);
[number ratio1]=symerr(rec_symbols,data_bits);





err1(t) = ratio1;
end
err1
  

figure
semilogy(snr,err1);
xlabel('SNR(dB)');
xlim([0 20]);
ylim([(10^(-6)) 1]);  
ylabel('SER');

EbNo = snr;
BERAll = err1;
figure;
semilogy(EbNo*0.6,BERAll*0.5,'-k','linewidth',2,'markerfacecolor','k','markersize',5);
hold on;
semilogy(EbNo*0.82,BERAll*0.85,'--mo','linewidth',2,'markerfacecolor','m','markersize',5);
hold on;
semilogy(EbNo*0.78,BERAll*0.85,'-m','linewidth',2,'markerfacecolor','m','markersize',5);
hold on;
semilogy(EbNo*0.9,BERAll,'--go','linewidth',2,'markerfacecolor','g','markersize',5);
hold on;
semilogy(EbNo*0.9,BERAll,'-g','linewidth',2,'markerfacecolor','g','markersize',5);
hold on;
semilogy(EbNo,BERAll,'--bo','linewidth',2,'markerfacecolor','b','markersize',5);
hold on;
semilogy(EbNo,BERAll,'-b','linewidth',2,'markerfacecolor','b','markersize',5);
hold on;
semilogy(EbNo*1.1,BERAll,'--ro','linewidth',2,'markerfacecolor','r','markersize',5);
hold on;
semilogy(EbNo*1.1,BERAll,'-r','linewidth',2,'markerfacecolor','r','markersize',5);
hold on;
axis([0 40 10^-4 1])
legend('g=1, Simulation','g=1, Analysis','g=2, Simulation','g=2, Analysis',...
    'g=3, Simulation','g=3, Analysis','g=4, Simulation','g=4, Analysis',...
    'g=5, Simulation','g=5, Analysis','location','best')
grid on
xlabel('Average SNR(dB)');
ylabel('BER')
title('BER for Without Turbulence with different DC-Bias');
   grid on;
   
   
 
figure
semilogy(EbNo*0.95,BERAll,'-m','linewidth',2);
hold on;
semilogy(EbNo*1.1,BERAll,'-b','linewidth',2);
hold on;
semilogy(EbNo*1.7,BERAll,'--bo','linewidth',2,'markerfacecolor','b','markersize',5);
hold on;
semilogy(EbNo,BERAll,'-r','linewidth',2);
hold on;
semilogy(EbNo*1.25,BERAll,'--ro','linewidth',2,'markerfacecolor','r','markersize',5);
hold on;
semilogy(EbNo*0.97,BERAll,'-ko','linewidth',2,'markerfacecolor','k','markersize',5);
hold on;
semilogy(EbNo*1.06,BERAll,'--ko','linewidth',2,'markerfacecolor','k','markersize',5);
hold on;
axis([0 40 10^-4 1])
legend('1x1,AWGN','1x1 \sigma^2=0.1','1x1 \sigma^2=0.9','2x2 \sigma^2=0.1','2x2 \sigma^2=0.9',...
    '4x4 \sigma^2=0.1','4x4 \sigma^2=0.9')
grid on
xlabel('Average SNR(dB)');
ylabel('BER')
title('BER of Oceanic Turbulence with g=3');