function [ thetaPh, phiPh, directionPh_new ] = scattering(g, directionPh)
% Scattering function 
% Input parameters: 
% g: HG asymmetry parameter 
% directionPh: photon’s current direction % 
% Returned parameters: 
% thetaPh: photon’s zenith angle after scattering 
% phiPh: photon’s azimuth angle after scattering 
% directionPh_new: photon’s direction after scattering 
% new azimuthal angle

phiPh = 2 * pi * rand; % new zenith angle 
x = rand; 
if (g == 0) 
    thetaPh = acos(2 * x - 1); 
else
    thetaPh = acos((1 + g^2 - ((1-g^2)/(1-g+2*g*x))^2) / 2*g); 
end
% new direction 
cos_phi = cos(phiPh); 
sin_phi = sin(phiPh); 
cos_theta = cos(thetaPh); 
sin_theta = sin(thetaPh); 
miX = directionPh(1);
miY = directionPh(2); 
miZ = directionPh(3); 
if (abs(miZ) > 0.99999) 
    miX_new = sin_theta * cos_phi; 
    miY_new = sin_theta * sin_phi; 
    miZ_new = sign(miZ) * cos_theta; 
else
    miX_new = (sin_theta / sqrt(1-miZ^2)) * (miX*miZ*cos_phi - miY*sin_phi) + miX*cos_theta; 
    miY_new = (sin_theta / sqrt(1-miZ^2)) * (miY*miZ*cos_phi + miX*sin_phi) + miY*cos_theta; 
    miZ_new = -sin_theta * sqrt(1-miZ^2) * cos_phi + miZ*cos_theta; 
end
directionPh_new = [miX_new, miY_new, miZ_new]; 
end