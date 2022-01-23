% Answer 2. )

% Group Delay of Second ORder Bandpass Filter using Butterworth Filter.
[z,p,k] = butter(6,0.2);
sos = zp2sos(z,p,k);


grpdelay(sos,128)

gd = grpdelay(sos,512);

[h,w] = freqz(sos,512);
pd = -unwrap(angle(h))./w;

plot(w/pi,gd,w/pi,pd)
grid
xlabel 'Normalized Frequency (\times\pi rad/sample)'
ylabel 'Group and phase delays'
legend('Group delay','Phase delay')




% Group Delay Response of Filter.
d = designfilt('lowpassiir','FilterOrder',6, ...
    'HalfPowerFrequency',0.2,'DesignMethod','butter');
grpdelay(d)


% Magnitude Response.
fs = 200;
d = designfilt('arbmagfir', ...
       'FilterOrder',88, ...
       'NumBands',4, ...
       'BandFrequencies1',[0 20], ...
       'BandFrequencies2',[25 40], ...
       'BandFrequencies3',[45 65], ...
       'BandFrequencies4',[70 100], ...
       'BandAmplitudes1',[2 2], ...
       'BandAmplitudes2',[0 0], ...
       'BandAmplitudes3',[1 1], ...
       'BandAmplitudes4',[0 0], ...
       'SampleRate',fs);
freqz(d,10:1/fs:78,fs)
filtord(d)

grpdelay(d,10:1/fs:78,fs)