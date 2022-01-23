function [ thetaPh, phiPh, positionPh, directionPh, weightPh ] = initialize_photon( thetaTr, phiTr, positionTr ) 
% Initialize the transmitted photon % Input parameters: 
% thetaTr: transmitter’s elevation angle 
% phiTr: transmitter’s aperture angle 
% positionTr: Transmitter’s position in (x, y, z) coordinates system 

% Returned parameters: 
% thetaPh: photon’s initial zenith angle 
% phiPh: photon’s initial azimuth angle 
% positionPh: Photon’s initial position in (x, y, z) coordinates system 
% directionPh: photon’s initial direction  
% weight: photon’s initial weight 
positionPh = positionTr; %Initial position of photon 
%phi of photon 
phiPh = 2*pi*rand; %random value in [0, 2?]
%theta of photon 
u = acos(1 - rand*(1 - cos(phiTr/2))); 
thetaPh = u + pi/2 - thetaTr; 
%direction of photon 
cos_phi = cos(phiPh); 
sin_phi = sin(phiPh); 
cos_theta = cos(thetaPh); 
sin_theta = sin(thetaPh); 
directionPh = [cos_phi*sin_theta, sin_phi*sin_theta, cos_theta]; 
%wheight of photon 
weightPh = 1e-5; 
end
