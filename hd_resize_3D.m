% -----------------------------------------------------------------------------
%  This function can be used for re-sizing 3D HD

% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Stanford University.
% ----------------------------------------------------------------------------*/

function [ HD2 ] = hd_resize_3D(HD, newsize)

HD2 = zeros(newsize(1), newsize(2), size(HD,3));

for L = 1:size(HD,3)
    hd_local = HD(:,:,L);
    HD2(:,:,L) = hd_resize_2D(hd_local,newsize);
end;

end

