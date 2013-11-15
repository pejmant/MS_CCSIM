% -----------------------------------------------------------------------------
% This function can be used for re-sizing 3D TI.

% Different resizing techniques can be used (e.g. Bicubic interpolation, 
% Nearest-neighbor interpolation, Bilinear interpolation)

% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Stanford University.
% ----------------------------------------------------------------------------*/

function [ ti2 ] = ti_resize_3D(ti, newsize)

ti2 = zeros(newsize(1), newsize(2), size(ti,3));

for L = 1:size(ti,3)
    ti_local = ti(:,:,L);
    ti2(:,:,L) = imresize(ti_local, [newsize(1) newsize(2)], 'nearest');
end;

end

