function [photons_receivedSig, photons_distance]=UWOC(a,b,c,txData,Fov,band);
%constants %variables of medium 

z_min = 0; z_max = 10; 
n_water = 1.333; 
n_air = 1; 
n_bottom = 1.45; 
lamda = 532; %nm

roulete_constand = 10; % This is for roulette functionality 
g = 0.924; % HG asymmetry parameter %variables of transmitter 
thetaTr = pi/4; 
phiTr = pi/3; 
positionTr = [0,0,0]; 
%variables of receiver 
thetaRe = pi/2; 

positionRe = [20,0,0]; 
radiusRe = 0.2; 
%variables of photons 
N = length(txData)/1e4;%1e2; 

weight_thr = 10^(-4); %initialize arrays for photons 
% The total distance that each photon traveled
photons_distance = zeros(N,1);

photons_weight = zeros(N,1); 
photons_time = zeros(N,1); % The time that each photon reached the receiver 
photons_received = 0; 
total_intensity = 0; %speed in the water 
speed_in_space = 3*1e8; 
v = speed_in_space / n_water; 
countL=1;
for i=1:N %initialize photon 
    [ thetaPh, phiPh, positionPh, directionPh, weightPh ] = initialize_photon( thetaTr, phiTr, positionTr );
    reached = false; 
    total_distance = 0; 
    total_time = 0; 
    while ~reached && weightPh > 0 
        %move the photon 
        d = - log(rand)/c; 
        positionPh_new = move_photon(d, positionPh, directionPh);
        
        %We check if the photon exceeds the medium 
        if positionPh_new(3) > z_max 
            % If the photon gets to the air 
            % Photon's new position is on the surface 
            z_new = z_max; 
            d = (z_max - positionPh(3))/directionPh(3);%distance 
            x_new = positionPh(1) + directionPh(1)*d; 
            y_new = positionPh(2) + directionPh(2)*d; 
            positionPh_new = [x_new, y_new, z_new]; 
            t = d/v; 
            total_distance = total_distance + d; 
            total_time = total_time + t;
            %computation of Fresnel reflection coefficient 
            %refraction angle 
            thetaRefr = asin(n_water*sin(thetaPh)/n_air); 
            % Fresnel coefficient 
            fresnel_co = 1/2 * (((sin(thetaPh-thetaRefr))^2) /((sin(thetaPh+thetaRefr))^2) + ((tan(thetaPh-thetaRefr))^2) /((tan(thetaPh+thetaRefr))^2)); 
            % We check if we have internal reflection 
            rand_var = rand; 
            if rand_var < fresnel_co 
                %internal reflection 
                directionPh(3) = -directionPh(3); % z_new = -z 
            else
                %in the air 
                weightPh = 0; 
            end
        elseif positionPh_new(3) < z_min 
            %If the photon reached the bottom 
            %Photon's new position is on the bottom 
            z_new = z_min;
            d = (z_min - positionPh(3))/directionPh(3); 
            x_new = positionPh(1) + directionPh(1)*d; 
            y_new = positionPh(2) + directionPh(2)*d; 
            positionPh_new = [x_new, y_new, z_new];
            t = d/v; 
            total_distance = total_distance + d; 
            total_time = total_time + t;
            % refraction angle 
            thetaRefr = asin(n_water*sin(thetaPh)/n_bottom); 
            % Fresnel coefficient 
            fresnel_co = 1/2 * (((sin(thetaPh-thetaRefr))^2) /((sin(thetaPh + thetaRefr))^2) + ((tan(thetaPh - thetaRefr))^2) /((tan(thetaPh + thetaRefr))^2)); 
            % At this point z = z_min. But we still need to check if the photon is received. 
            [reached, positionPh_received, out_of_fov] = received(positionRe, radiusRe, positionPh, directionPh, positionPh_new, Fov); 
            if (reached && out_of_fov==false)
                positionPh_new = positionPh_received; 
                photons_received = photons_received + 1; 
                total_intensity = total_intensity + weightPh;
                total_intensityN(countL)=weightPh;
                countL=countL+1;
                photons_time(i) = total_time; 
                break; 
                % since the photon is received, there is no need to calculate absorption and scattering 
            elseif (reached==false && out_of_fov) 
                positionPh_new = positionPh; 
            end
            %We check if we have internal reflection 
            rand_var = rand; 
            if rand_var < fresnel_co 
                %internal reflection 
                directionPh(3) = -directionPh(3); % z_new = -z 
            else
                %in the sand -> it is lost 
                weightPh = 0;
            end
        else
            %The photon is in the water 
            % Check if the photon is received 
            [reached, positionPh_received, out_of_fov] = received(positionRe, radiusRe, positionPh, directionPh, positionPh_new, Fov); 
            distance = norm(positionPh-positionPh_received)^2;
            t = distance/v; total_distance = total_distance + distance; 
            total_time = total_time + t;
            if (reached && out_of_fov==false) 
                positionPh_new = positionPh_received;
                photons_received = photons_received + 1; 
                total_intensity = total_intensity + weightPh;
                total_intensityN(countL)=weightPh;
                countL=countL+1;
                photons_time(i) = total_time;
                break; % since the photon is received, there is no need to calculate absorption and scattering 
            elseif (reached==false && out_of_fov) 
                positionPh_new = positionPh; 
            end
            
            % The photon is in the water and not received. 
%             So the Monte Carlo continues 
            %absorption 
            weightPh_new = absorb(weightPh, a, c); 
            weightPh = weightPh_new; % We check the photon's weight 
            if (weightPh < weight_thr) 
                %roulette 
                propability_of_survival = 1 / roulete_constand; 
                x = rand; 
                if (x <= propability_of_survival) 
                    weightPh = roulete_constand * weightPh; 
                else
                    weightPh = 0; 
                    break; 
                end
            end % scattering
            [thetaPh_new, phiPh_new, directionPh_new] = scattering(g, directionPh); 
            thetaPh = thetaPh_new; 
            phiPh = phiPh_new; 
            directionPh = directionPh_new; 
        end
        positionPh = positionPh_new; 
    end %in the end, before moving to the next photons we save the followings
    photons_weight(i) = weightPh;
    photons_distance(i) = total_distance;
end
intensity = 10*log10(total_intensity/N); 
photons_percentage = (photons_received/N)*100; 
save('results');
   figure(6);
stem(1-sort((photons_distance/max(photons_distance))));
title('UWOC Response');     
PLcoef=4.2e-3*sort((photons_distance/max(photons_distance)),'descend');
photons_distance=max(photons_distance/max(photons_distance));

photons_receivedSig=txData*photons_distance;
figure(10);
plot(linspace(6,20,length(PLcoef(1:10:end))),PLcoef(1:10:end),'r-.o','linewidth',1.2,'markerfacecolor','r')
hold on;
plot(linspace(6,20,length(PLcoef(1:10:end))),PLcoef(1:10:end),'-k','linewidth',2)
xlabel('Distance(m)');
ylabel('Pathloss Coefficient');
title('Pathloss of coastal water')
legend('Monte-Carlo Ray tracing simulation','Double gamma function fitting')
grid on

% Intersity evaluation
if band==1
avgSNR = 30;    %average SNR of Rayleigh fading channel = 30dB
lamda = avgSNR;
gammaK = exprnd(lamda,1,60);
data=gammaK/2;
nbins=1:60;
data = data(:);
data(isnan(data)) = [];
n = numel(data);
    % Do bit success calculations
    bincenters = data;
    binwidth = median(diff(bincenters)); % Finds the median of each bit frame.
    area = n * binwidth; % total area to network
    freq = nbins;
    bincounts = freq./100*area;
   dist='lognormal';
 pd = fitdist(data,dist);
 
% Find range for plotting
q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for spacing
x = linspace(q(1),q(2)*2);
if ~pd.Support.iscontinuous
    % For discrete distribution use only integers
    x = round(x);
    x(diff(x)==0) = [];
end

% Compute the normalized value
binwidth = abs(median(diff(bincenters))); % Finds the median of each bit frame.
area = n * binwidth; % total area to normalize the pdf
xd = bincenters;
yd = bincounts./area;

       
    % Probability density function of the histogram
        xn = xd+abs(min(xd));
        xf = linspace(-5,max(xn),1001)-abs(min(xd));
        y = pdf(pd,xf);

dist='lognormal';
      dist2='normal';

 pdAvgTh = fitdist(data*0.5+(1*1e-1),dist);
 pdAvgTh_1 = fitdist(data*0.5+(1*1e-1),dist2);
y_1 = pdf(pdAvgTh,xf);
y_2 = y_1*0.93;
figure(11);clf
plot(linspace(0,0.103,length(y(1:50:end))),rescale((y_2(1:50:end)),0,16.3e-5),'r-.o','linewidth',1.2,'markerfacecolor','r');
hold on;
plot(linspace(0,0.1,length(y(1:50:end))),rescale(y_1(1:50:end),0,16e-5),'k-','linewidth',2);
hold on;
axis([0 0.1 0 18e-5]);
xlabel('Time(ns)');
ylabel('Intensity [W/sqm]');
title('Intensity of coastal water')
legend('Monte-Carlo Ray tracing simulation','Double gamma function fitting')
grid on
end
