
% -----------------------------------------------------------------------------
% This function can be used for both 3D unconditional and conditional
% MS-CCSIM. For using in unconditional mode, the hd should be just NaN
% value.

% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for Reservoir Forecasting, Stanford University.
% ----------------------------------------------------------------------------*/


function [ C, error_location ] = CCSIM_3D(ti, hd, T, OL, CT, fc, prop, rad, cand, mrp, T_vibration, ms_level, real_numb)


%% Input Parameters
% - ti: Training image
% - hd: Hard data matrix
% - T: Template size
% - OL: Overlap size
% - CT: Size of Co-Template (e.g. [2 2])
% - fc: Number of facies (if fc = 0, histogram matching will not apply; default is 0)
% - prop: The proportion of TI which will be scanned to find the matched hard data (0 < prop <=1))
% - rad: Neighbourhood radius for MS simulation (default is [1 1])
% - cand: Number of candidates for pattern pool
% - mrp: Multiple random path flag (1:on, 0:off).
% - T_vibration : Multiple template size flag (1:on, 0:off).
% - ms_level: Number of MS level (1, 2, 3).
% - real_numb: Number of the realizations

%% Output Parameters
% - C: Simulation grid for output
% - error_location: The locations where the HD is macthed (all zero means
% no mismatch)

%%
if ms_level >=1
    HD1_1 = hd; ti1_1 = ti; lag = 1;
    if mrp == 1
        lag = 6;
        HD2_1 = flipdim(HD1_1,1); ti2_1 = flipdim(ti1_1,1);
        HD3_1 = flipdim(HD1_1,2); ti3_1 = flipdim(ti1_1,2);
        HD4_1 = flipdim(flipdim(HD1_1,1),2); ti4_1 = flipdim(flipdim(ti1_1,1),2);
        HD5_1 = transpose_3D(HD1_1); ti5_1 = transpose_3D(ti1_1);
        HD6_1 = transpose_3D(HD3_1); ti6_1 = transpose_3D(ti3_1);
    end;
end;

C = zeros(numel(hd),1);

if ms_level >=2
    HD1_2 = hd_resize_3D(HD1_1, size(HD1_1)/2); ti1_2 = ti_resize_3D(ti1_1, size(ti1_1)/2);
    if mrp == 1
        HD2_2 = hd_resize_3D(HD2_1, size(HD2_1)/2); ti2_2 = ti_resize_3D(ti2_1, size(ti2_1)/2); 
        HD3_2 = hd_resize_3D(HD3_1, size(HD3_1)/2); ti3_2 = ti_resize_3D(ti3_1, size(ti3_1)/2); 
        HD4_2 = hd_resize_3D(HD4_1, size(HD4_1)/2); ti4_2 = ti_resize_3D(ti4_1, size(ti4_1)/2);
        HD5_2 = hd_resize_3D(HD5_1, size(HD5_1)/2); ti5_2 = ti_resize_3D(ti5_1, size(ti5_1)/2); 
        HD6_2 = hd_resize_3D(HD6_1, size(HD6_1)/2); ti6_2 = ti_resize_3D(ti6_1, size(ti6_1)/2);
    end;
end;
    
    
if ms_level >=3
    HD1_3 = hd_resize_3D(HD1_1, size(HD1_1)/4); ti1_3 = ti_resize_3D(ti1_1, size(ti1_1)/4);
    if mrp == 1
        HD2_3 = hd_resize_3D(HD2_1, size(HD2_1)/4); ti2_3 = ti_resize_3D(ti2_1, size(ti2_1)/4); 
        HD3_3 = hd_resize_3D(HD3_1, size(HD3_1)/4); ti3_3 = ti_resize_3D(ti3_1, size(ti3_1)/4); 
        HD4_3 = hd_resize_3D(HD4_1, size(HD4_1)/4); ti4_3 = ti_resize_3D(ti4_1, size(ti4_1)/4);
        HD5_3 = hd_resize_3D(HD5_1, size(HD5_1)/4); ti5_3 = ti_resize_3D(ti5_1, size(ti5_1)/4);
        HD6_3 = hd_resize_3D(HD6_1, size(HD6_1)/4); ti6_3 = ti_resize_3D(ti6_1, size(ti6_1)/4);
    end;
end;
    
    
% % h = waitbar(0, sprintf('CCSIM is running ... %i realization(s)', real_numb), ...
% %     'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
% % setappdata(h,'canceling',0)


R1 = 1:lag:real_numb;
R2 = 0; R3 = 0; R4 = 0; R5 = 0; R6 = 0;
if mrp == 1
    R2 = 2:lag:real_numb;
    R3 = 3:lag:real_numb;
    R4 = 4:lag:real_numb;
    R5 = 5:lag:real_numb;
    R6 = 6:lag:real_numb;
end


if (real_numb > 1) && (T_vibration == 1)
    T_new = ones(real_numb,3);
    T_new(1:ceil(real_numb/3),1) = T(1);
    T_new(1:ceil(real_numb/3),2) = T(2);
    T_new(1:ceil(real_numb/3),3) = T(3);
    
    T_new(1+ceil(real_numb/3):2*ceil(real_numb/3),1) = T(1) + 4;
    T_new(1+ceil(real_numb/3):2*ceil(real_numb/3),2) = T(2) + 4;
    T_new(1+ceil(real_numb/3):2*ceil(real_numb/3),3) = T(3);    
    
    T_new(1+2*ceil(real_numb/3):end,1) = T(1) - 4;
    T_new(1+2*ceil(real_numb/3):end,2) = T(2) - 4;
    T_new(1+2*ceil(real_numb/3):end,3) = T(3);    
    
    T_new(:,4) = rand(size(T_new,1),1);
    T_new = sortrows(T_new,4);
    T_new = T_new(:,1:3);
else
    T_new = repmat(T,[real_numb 1]);
end;


error_location = zeros(numel(find(~isnan(hd))),real_numb);

for i = 1:real_numb
    
%     if getappdata(h,'canceling')
%         break
%     end     
    
    if any(R1==i)
        hd0 = HD1_1; ti0 = ti1_1;
        if ms_level >=2
            hd1 = HD1_2; ti1 = ti1_2;
        end;
        if ms_level >=3
            hd2 = HD1_3; ti2 = ti1_3;
        end;
        
    elseif any(R2==i)
        hd0 = HD2_1; ti0 = ti2_1;
        if ms_level >=2
            hd1 = HD2_2; ti1 = ti2_2;
        end;
        if ms_level >=3
            hd2 = HD2_3; ti2 = ti2_3;        
        end;
        
    elseif any(R3==i)
        hd0 = HD3_1; ti0 = ti3_1;
        if ms_level >=2
            hd1 = HD3_2; ti1 = ti3_2;
        end;
        if ms_level >=3
            hd2 = HD3_3; ti2 = ti3_3;
        end;
        
    elseif any(R4==i)
        hd0 = HD4_1; ti0 = ti4_1;
        if ms_level >=2
            hd1 = HD4_2; ti1 = ti4_2;
        end;
        if ms_level >=3
            hd2 = HD4_3; ti2 = ti4_3;
        end;
        
    elseif any(R5==i)
        hd0 = HD5_1; ti0 = ti5_1;
        if ms_level >=2
            hd1 = HD5_2; ti1 = ti5_2;
        end;
        if ms_level >=3
            hd2 = HD5_3; ti2 = ti5_3;
        end;
        
    elseif any(R6==i)
        hd0 = HD6_1; ti0 = ti6_1;
        if ms_level >=2
            hd1 = HD6_2; ti1 = ti6_2;
        end;
        if ms_level >=3
            hd2 = HD6_3; ti2 = ti6_3;        
        end;
    end;    
    

    tStart = tic;
    if ms_level ==1
        [MS0, ~] = CCSIM_3D_MS2(ti0, hd0, T_new(i,:), OL, CT, fc, prop, cand);
    elseif ms_level==2
        [~, LOC1] = CCSIM_3D_MS2(ti1, hd1, [T_new(i,1)/2 T_new(i,2)/2 T(3)],...
            [OL(1)/2 OL(2)/2 OL(3)], CT, fc, prop, cand);
        [MS0, ~] = CCSIM_3D_MS1(ti0, hd0, LOC1, T_new(i,:), OL, fc, rad);
    else
        [~, LOC2] = CCSIM_3D_MS2(ti2, hd2, [T_new(1)/4 T_new(2)/4 T(3)], [OL(1)/4 OL(2)/4 OL(3)], CT, fc, prop, cand);
        [~, LOC1] = CCSIM_3D_MS1(ti1, hd1, LOC2, [T_new(i,1)/2 T_new(i,2)/2 T(3)], [OL(1)/2 OL(2)/2 OL(3)], fc, rad);
        [MS0, ~] = CCSIM_3D_MS1(ti0, hd0, LOC1, T_new(i,:), OL, fc, rad);
    end;
    
    tEnd = toc(tStart);
    
    
    
    if any(R2==i)
        MS0 = flipdim(MS0,1); 
                
    elseif any(R3==i)
        MS0 = flipdim(MS0,2);
                
    elseif any(R4==i)
        MS0 = flipdim(flipdim(MS0,1),2);
        
    elseif any(R5==i)     
        MS0 = transpose_3D(MS0);
        
    elseif any(R6==i)
        MS0 = transpose_3D(MS0);
    else
        MS0 = MS0(1:size(hd,1),1:size(hd,2),1:size(hd,3));          
    end;
    
  
    
    C(:,i) = MS0(:);
    
    hd_all = numel(find(~isnan(hd0)));
    
    if hd_all ~=0
        [error, error_location(:,real_numb)] = hd_error(hd,MS0);
        mis_hd = 100*(error/hd_all);
        disp(['** Mismatch HD: ',num2str(mis_hd),'% ***'])
    end;
    

    disp(['********  CPU time for the grid size of ',num2str(size(hd,1)),'x',...
        num2str(size(hd,2)), 'x',num2str(size(hd,3)), ' is ', num2str(tEnd),...
        '  ********'])
    
    
    
% waitbar(i / real_numb, h, sprintf('CCSIM is running...Please wait...%ith realization completed', i))
    
end;

% delete(h) 

end

