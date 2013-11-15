
% -----------------------------------------------------------------------------
% This function can be used for finding the best minimum error boundary

% Authors: Efros, A. A., and W. T. Freeman
% Reference: Image quilting for texture synthesis and transfer paper 
% presented at Proceedings of the ACM SIGGRAPH 
% Conference on Computer Graphics

%function Y = mincut(X,dir)
%Computes the minimum cut from one edge of the matrix to the other
%
%Inputs:
%   X:  Evalutations of 2d function to be cut along local minima or overlap
%   region
%   dir: 0 = vertical cut, 1 = horizontal cut
%
%Outputs:
%   C: Matrix containing entries indicating side
%       -1: left (or top) side of cut
%        0: along the cut
%        1: right (or bottom) side of cut
% ----------------------------------------------------------------------------*/

function C = mincut(X,dir)

%-------------------------------------------------------------------------------

if( nargin > 1 && dir == 1 )
    X = X';
end;

%---------------------------------------------------------------------------------
%Allocate the current cost array, and set first row to first row of X

E = zeros(size(X));
E(1:end,:) = X(1:end,:);

%Starting with the second array, compute the path costs until the end
for i=2:size(E,1),
    
    E(i,1) = X(i,1) + min( E(i-1,1), E(i-1,2) );
    for j=2:size(E,2)-1,
        E(i,j) = X(i,j) + min( [E(i-1,j-1), E(i-1,j), E(i-1,j+1)] );
    end;
    E(i,end) = X(i,end) + min( E(i-1,end-1), E(i-1,end) );
end;

%---------------------------------------------------------------------------------------------
%Backtrace to find the cut

C = zeros(size(X));
[cost, idx] = min(E(end, 1:end));
C(i, 1:idx-1) = -1;
C(i, idx) = 0;
C(i, idx+1:end) = +1;  
%----------------------------------------------------------------------------------------------       
for i=size(E,1)-1:-1:1,
    for j=1:size(E,2),

        if( idx > 1 && E(i,idx-1) == min(E(i,idx-1:min(idx+1,size(E,2))) ) )
            idx = idx-1;
        elseif( idx < size(E,2) && E(i,idx+1) == min(E(i,max(idx-1,1):idx+1)) )
            idx = idx+1;
        end;
               
        C(i, 1:idx-1) = -1;
        C(i, idx) = 0;
        C(i, idx+1:end) = +1;            
    end;
end;
%--------------------------------------------------------------------------------------------------   
if( nargin > 1 && dir == 1 )
    C = C';
end;



