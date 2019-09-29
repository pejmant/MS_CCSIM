function [ out ] = transpose_3D( in )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


out = zeros(size(in,2),size(in,1),size(in,3),size(in,4));
for i = 1:size(in,4)
    for j = 1:size(in,3)
        a = in(:,:,j,i);
        a = a';
        out(:,:,j,i) = a;
    end;
end
