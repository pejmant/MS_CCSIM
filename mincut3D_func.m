
% -----------------------------------------------------------------------------
%  This function can be used for finding the best minimum error boundary in
%  3D simulations.

% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Stanford University.
% ----------------------------------------------------------------------------*/

function [ M ] = mincut3D_func( target_final, imout, T, OL, i, j, k )

M = ones(T(1),T(2),T(3));

if (j>1)
    E = ( target_final(1:1+T(1)-1,1:1+OL(2)-1,1:1+T(3)-1) - imout(1:1+T(1)-1,1:1+OL(2)-1,1:1+T(3)-1) ).^2;
    
    for L = 1:T(3)
        % Left OL 
        E_now = E(:,:,L);
        % Compute the mincut array
        C = mincut(E_now, 0);
        % Compute the mask and write to the destination
        M(1:end, 1:OL(2), L) = double(C >= 0);
    end
end


if (i>1)
    E = ( target_final(1:1+OL(1)-1, 1:1+T(2)-1, 1:1+T(3)-1) - imout(1:1+OL(1)-1, 1:1+T(2)-1, 1:1+T(3)-1) ).^2;
    
    for L = 1:T(3)
        % We have a bottom OL
        E_now = E(:,:,L);
        % Compute the mincut array
        C = mincut(E_now, 1);
        % Compute the mask and write to the destination
        M(1:OL(1), 1:end, L) = M(1:OL(1), 1:end, L) .* double(C >= 0);
    end;
end

if (k>1)
    E = ( target_final(1:1+T(1)-1, 1:1+T(2)-1, 1:1+OL(3)-1) - imout(1:1+T(1)-1, 1:1+T(2)-1, 1:1+OL(3)-1) ).^2;
    
    for L = 1:T(2)
        % Left OL 
        E_now = reshape(E(:,L,:),size(E,1),size(E,3));
        % Compute the mincut array
        C = mincut(E_now, 0);
        C = reshape(C,size(E,1),1,size(E,3));
        % Compute the mask and write to the destination
        M(1:T(1), L, 1:OL(3)) = M(1:T(1), L, 1:OL(3)) .*double(C >= 0);
    end
end
end