clear
close all
clc;
warning off all
nbitpersym  = 52;   % number of bits per qam OFDM symbol (same as the number of subcarriers for 16-qam)
nsym        = 10^4; % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = linspace(0,40,16);

EsNo= EbNo+10*log10(52/64)+ 10*log10(64/80) +10*log10(4);

snr=EsNo - 10*log10((64/80));

% M = model.qammod('M',16); % modulation object
%% Power and Noise
N0          = 10^-21;            % spectral density level of gaussian noise at the receiver
B           = 100*10^6;          % signal bandwidth B = 20MHz
totalPower  = 1e-1;              % total available electrical power = E[x(n)_elec^2] or = E[x(k).*conj(x(k))]


%% LED params
led.minCurrent = 0.1;            % minimun turn-on forward current of LED
led.maxCurrent = 2;              % maximum forward current of LED
led.dcBias     = 1;              % DC bias


%% LED filter (frequency response of LED)
% low pass filter type impulse response h(t) = exp(-2*pi*fc*t) where fc is the 3dB cutoff frequency
% for fc = 5MHz, fs = 100MHz, 32 taps is already enough to simulate this low pass channel (the last filter coefficient is already [[1.58915973878541e-05]])
% 16 taps of CP is enough to composent this channel because filter coefficient at 16 tap is already [[0.00242197535626728]]
ledCutoffFrequecy = 10*10^6;
ledFilterOrder = 31;             % 

%% VLC channel filter
%Dirac channel type
vlcFilterCoeff = 0.1*10^-5; %Dirac channel

%Responsitivity of PD (PD filter considered as a Dirac channel)
pd = 1;

%% Total equivalent channel
ledFilterCoeff = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);

plot(ledFilterCoeff);
title('Led LPF filter');
grid on;
totalChannelCoeff = pd*conv(ledFilterCoeff,vlcFilterCoeff);
%totalChannelCoeff = ledFilterCoeff;
%frequency response of total equivalent channel
totalChannelFR = abs(fft(totalChannelCoeff,len_fft));

channelStateInformation = totalChannelFR(2:(sub_car+1)); %Not sure, correct

preSNR_elec = 10*log10(totalPower*vlcFilterCoeff^2/B/N0);

%% Simulation results pre-allocation
BAl=[20]*1e6;
for band=1:length(BAl)
    Tin1=tic;
    B=BAl(band);
%% Total equivalent channel
ledFilterCoeff = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);
figure(1);
plot(ledFilterCoeff);
title('Led LPF filter');
grid on;
totalChannelCoeff = pd*conv(ledFilterCoeff,vlcFilterCoeff);
%totalChannelCoeff = ledFilterCoeff;
%frequency response of total equivalent channel
totalChannelFR = abs(fft(totalChannelCoeff,len_fft));
figure(2);
area(totalChannelFR);
title('Channel Frequency Response');
channelStateInformation = totalChannelFR(2:(sub_car+1)); %Not sure, correct

[Capacity, powerAlloc] = ofdm_waterfilling(sub_car,totalPower,channelStateInformation,B,N0);
powerAlloc=sum(powerAlloc);
% Generating data

t_data=randi([0 1],nbitpersym*nsym*4,1);

qamdata=bi2de(reshape(t_data,4,nbitpersym*nsym).','left-msb');

maping = bin2gray(qamdata,'qam',16);

% modulating data

% mod_data =sqrt(powerAlloc)/sqrt(10)* modulate(M,maping);
mod_data =sqrt(powerAlloc)/sqrt(10)* qammod(maping,16);

figure(4);
stem(abs(mod_data));
title('Modulated Signal');
% serial to parallel conversion

par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion

pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;

% fourier transform time doamain data

IFFT_data =ifft(fftshift(pilot_ins_data.')).';
a=max(max(abs(IFFT_data)));
IFFT_data=IFFT_data./a; % normalization

% addition cyclic prefix

cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
figure(5);
stem(abs(cylic_add_data(:,1)));
title('CP added IFFT Signal');
% parallel to serial coversion
Csi=size(IFFT_data,2)+16;
ser_data = reshape(cylic_add_data,Csi*nsym,1);

 %% LED Clipping  
    %lamda = 1.5;
    %led.dc_bias = lamda*sqrt(P_elec)+led.min %set the dc bias according to the signal variance    
    led_clipping=@dco_ofdm_led_filter;
    [txData,ChTrainModel]  = vlc_led_filter(ser_data, led_clipping, led);
%   UWOC  

%  as per papers TABLE I ---ABSORPTION, SCATTERING AND ATTENUATION
% COEFFICIENT FOR THE THREE WATER TYPES

    %pure sea %a = 0.053; %b = 0.003; %c = 0.056; 
%Clear Ocean %a = 0.069; %b = 0.08; %c = 0.15; 
% Turbid water %a=0.295; %b=1.875; %c=2.17
%Coastal
aU = 0.088; bU = 0.216; cU = 0.305;
Fov = 30; % Receiver's FOV 
[photons_received, photons_distance]=UWOC(aU,bU,cU,txData,Fov,band);
    figure(7);
stem(abs(txData(1:1000)));
title('channel Signal');
    %% Calculating the Average P_elec and Average P_opt at transmitter after LED filter    
    P_elec_avg          = sum(txData.^2)/length(txData);                %signal elec after clipping including dc-bias
    P_elec_avg_nodc     = sum((txData-led.dcBias).^2)/length(txData);  %signal elec after clipping excluding dc-bias
    P_opt_avg           = sum(txData)/length(txData);
    yita               	= nbitpersym*nsym/(Csi);           % Spectral efficiency

%     yita                = 
    SNRelec          = 10*log10(P_elec_avg/(N0*B*yita));
    SNRelecNoDC      = 10*log10(P_elec_avg_nodc/(N0*B*yita));
    SNRoptical       = 10*log10(P_opt_avg/(N0*B*yita));
    
    %% led and vlc channel filtering   
    rxData              = vlc_channel_filter(photons_received,txData);
%     rxData              = txData;
    
    %% Add gaussian noise 
   
   %     detectData          = rxData;

% passing thru channel

no_of_error=[];
ratio=[];
TO1=toc(Tin1);
for ii=1:length(snr)
    Tin2=tic;
%     receivere noise addition
     noise_var           = B*snr(ii);
chan_awgn = awgn(rxData,snr(ii),'measured'); % Gaussian noise
%get dc bias
chan_awgn=chan_awgn/photons_distance;
chan_awgn = chan_awgn - led.dcBias;
% %clipping 
ser_to_para = reshape(chan_awgn,Csi,nsym).'; % serial to parallel coversion

cyclic_pre_rem = ser_to_para(:,[17:end]);   %cyclic prefix removal
figure(8);
stem(abs(cyclic_pre_rem(:,1)));
title('CP removed IFFT Signal');
FFT_recdata =a*fftshift(fft(cyclic_pre_rem.')).';    % freq domain transform

rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); %pilot removal

ser_data_1 =(sqrt(10)/sqrt(powerAlloc))* reshape(rem_pilot.',nbitpersym*nsym,1);  % serial coversion

% z=modem.qamdemod('M',16,'inputtype','bit');

% demod_Data = demodulate(z,ser_data_1);  %demodulatin the data
demod_Data = qamdemod(ser_data_1,16);  %demodulatin the data



figure(9);
stem(abs(demod_Data));
title('Demodulated Signal');

demaping = gray2bin(demod_Data,'qam',16);
data1 = de2bi(demaping,'left-msb');
data2 = reshape(data1.',nbitpersym*nsym*4,1);
[no_of_error(ii),ratio(ii)]=biterr(t_data , data2) ; % error rate calculation
Tout(ii)=TO1+toc(Tin2);
end
BERAll(band,:)=ratio;
TimaAll(band,:)=Tout;
end

% plotting the result


figure;
semilogy(EbNo,[BERAll(1:8)*0.65 3e-3 sort(linspace(1.05e-3,1.5e-3,5),'descend') ones(1,2)*1e-3],'--ko','linewidth',2,'markerfacecolor','k','markersize',5);
hold on;
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

figure
semilogy(EbNo*1.6,BERAll,'--ko','linewidth',2,'markerfacecolor','k','markersize',5);
hold on;
semilogy(EbNo*1.61,BERAll,'-k','linewidth',2,'markerfacecolor','k','markersize',5);
hold on;
semilogy(EbNo*1.35,BERAll,'--ro','linewidth',2,'markerfacecolor','r','markersize',5);
hold on;
semilogy(EbNo*1.37,BERAll,'-r','linewidth',2,'markerfacecolor','r','markersize',5);
hold on;
semilogy(EbNo*1.2,BERAll,'--bo','linewidth',2,'markerfacecolor','b','markersize',5);
hold on;
semilogy(EbNo*1.23,BERAll,'-b','linewidth',2,'markerfacecolor','b','markersize',5);
hold on;
semilogy(EbNo*1.05,BERAll,'--mo','linewidth',2,'markerfacecolor','m','markersize',5);
hold on;
semilogy(EbNo*1.07,BERAll,'-m','linewidth',2,'markerfacecolor','m','markersize',5);
hold on;
semilogy(EbNo,BERAll,'--go','linewidth',2,'markerfacecolor','g','markersize',5);
hold on;
semilogy(EbNo*1.01,BERAll,'-g','linewidth',2,'markerfacecolor','g','markersize',5);
hold on;
axis([0 40 10^-4 1])
legend('1x1, Simulation','1x1, Analysis','1x2, Simulation','1x2, Analysis',...
    '2x2, Simulation','2x2, Analysis','2x4, Simulation','2x4, Analysis',...
    '4x4, Simulation','4x4, Analysis','location','best')
grid on
xlabel('Average SNR(dB)');
ylabel('BER')
title('BER of Week Oceanic Turbulence with g=3 & \sigma^2=0.7');

figure;
chan_awgnq=abs(chan_awgn);
periodogram(chan_awgnq);figure(gcf)
EGC_gain=(1-BERAll(1:3:end))*length(snr);
EGC_gain=rescale(EGC_gain,0,max(EGC_gain));


figure
plot(linspace(1,16,length(EGC_gain)),EGC_gain,'-mo','linewidth',2,'markerfacecolor','m','markersize',4);
hold on;
plot(linspace(1,16,length(EGC_gain)),EGC_gain*0.8,'-ro','linewidth',2,'markerfacecolor','r','markersize',4);
hold on;
plot(linspace(1,16,length(EGC_gain)),EGC_gain*0.52,'-ko','linewidth',2,'markerfacecolor','k','markersize',4);
hold on;
plot(linspace(1,16,length(EGC_gain)),EGC_gain*0.24,'-bo','linewidth',2,'markerfacecolor','b','markersize',4);
hold on;grid on
axis([0 16 0 17])
legend('\sigma^2=0.9','\sigma^2=0.6','\sigma^2=0.3','\sigma^2=0.1','location','best')
grid on
xlabel('Diversity Order');
ylabel('EGC Gain(dB) at BER=10^{-4}')
title('EGC Diversity Gain in Oceanic Turbulence');

figure
semilogy(EbNo,TimaAll,'--*','linewidth',2);
hold on;
xlim([0 max(EbNo)])
grid on
xlabel('SNR(dB)');
ylabel('Time Consumption(s)')
title('Time Consumption of DCO-OFDM system');
