function [ error, location ] = hd_error(hd,rel)

f = find(isnan(hd)==0);
error = 0;
location = zeros(1);

for i = 1: numel(f)
    if (hd(f(i))~=rel(f(i)))
        error = error + 1;
        location (error,1) = f(i);
    end;    
end

