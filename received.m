function [ boolean_value, positionPh_received, out_of_fov] = received ( positionRe, radiusRe, positionPh, directionPh, positionPh_new, fov) 
% Returns true if the photon is considered as received. We check 
% the cases where the photon may be received % 
% Input parameters: 
% positionRe: receiver’s position in (x, y, z) 
% radiusRe: receiver’s radius 
% positionPh: photon’s previous position 
% directionPh: photon’s current direction
% positionPh_new: photon’s current position 
% fov: receiver’s FOV % Returned parameters: 
% Boolean_value: this is true if the photon is considered as received. Otherwise, this is false 
% positionPh_received: the position at which the photon is received by the receiver 
% out_of_fov: This is true if the photon is considered as received (boolean_value = true) 
% from a geometrical aspect, but it is out of the receiver’s FOV, 
% so it will not be considered as received as the receiver doesn’t “see” it.
boolean_value = false; 
positionPh_received = positionPh_new; 
out_of_fov = false; 
if (((positionPh_new(1)-positionRe(1))^2 + (positionPh_new(2)- positionRe(2))^2 +(positionPh_new(3)-positionRe(3))^2)<= radiusRe^2) 
    boolean_value = true; 
else
    a = norm(directionPh)^2;
    b= 2*dot(directionPh,positionPh-positionRe); 
    c= norm(positionPh-positionRe)^2 - radiusRe^2; 
    discr = b^2 - 4*a*c; 
    if (discr>=0) 
        % This means that line of the photon's direction intersects with the sphere. 
%         Now we need to check if the point of the photon's new position on the line is before or after the point 
%         where it intersects with the sphere.
% For this reason, first we find the distance from the photon's 
% initial position with point(s) where the line intersects with the medium 
d1 = (-b + sqrt(discr)) / 2*a; 
d2 = (-b - sqrt(discr)) / 2*a; 
d_min = min(d1, d2); 

% Distance of phton's new position to its old position 
distancePh_new = norm(positionPh-positionPh_new); 
if (distancePh_new >= d_min) 
    % The photon passes inside the sphere, so it is received 
    boolean_value = true; 
    % Now the photon's new position must change so that it would be inside the receiver's aperture area. 
%     The position at which the photon is received is at distance d_min from its initial point 
positionPh_received = positionPh+directionPh*d_min; 
end
    end
end
if boolean_value == true 
    % We check the diarection's angle and FOV 
    line1 = positionPh(3) - positionRe(3); % ?Z = Zp - Zr 
    line2 = norm(positionPh - positionRe); % ?? 
    theta = acos(line1/line2); 
    if theta <= fov/2 
        boolean_value = true; 
        out_of_fov = false; 
    else
        boolean_value = false; 
        out_of_fov = true; 
    end
end
end
