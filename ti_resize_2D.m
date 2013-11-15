% -----------------------------------------------------------------------------
% This function can be used for re-sizing 2D TI.

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

function [ ti2 ] = ti_resize_2D(ti, newsize)

ti2 = imresize(ti, newsize, 'nearest');

end

