function [ position_new ] = move_photon( d, position_old, direction ) 
% Moves the photon to its new position according to d, old position and direction 
% Input parameters:
% d: step size ?d
% position_old: photon’s previous position 
% direction: photon’s current direction 
% Returned parameters:
% position_new: Photon’s new position 
position_new = position_old + direction*d; 
end
