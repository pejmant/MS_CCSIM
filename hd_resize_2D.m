% -----------------------------------------------------------------------------
% This function can be used for re-sizing 2D HD

% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Stanford University.
% ----------------------------------------------------------------------------*/

function [ HD2 ] = hd_resize_2D(HD, newsize)

% Hard data resizing
[a,b] = size(HD);

% find the non-NaN entries in HD
idx = ~isnan(HD);

% find their corresponding row/column indices
[i,j] = find(idx);

% resize your matrix as desired, i.e. scale the row/column indices
i = ceil(i*newsize(1)/a);
j = ceil(j*newsize(2)/b);

% write the old non-NaN entries to HDnew using accumarray
% you have to set the correct size of Gnew explicitly
% maximum value is chosen if many entries share the same scaled i/j indices
% NaNs are used as the fill
HD2 = accumarray([i, j], HD(idx), [newsize(1) newsize(2)], @max, NaN);

end

