
% -------------------------------------------------------------------------
% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Energy Resources Engineering 
% Department,Stanford University.
% -----------------------------------------------------------------------*/



function [Grid_Sim, LOC] = hd_sp_3D_MS(TI, hd, T, OL, CT, fc, prop, cand, i, j, k, grid)


%% Input Parameters
% - TI: Training image
% - hd: hard data matrix
% - sizeout: Size of output that should be considered larger than its actuall size 
% - T: Template size
% - OL: Overlap size
% - cand: number of candidates
% - prop: Proportion of the TI that should be scanned for HD [0 1]


%% Output Parameters
% - Grid_Sim: Simulation grid for output

%% ---------------------------------------------------------------------------------------- 
sizeout_orig = size(hd);

sizeout(1) = sizeout_orig(1) + 2*(T(1));
sizeout(2) = sizeout_orig(2) + 2*(T(2));
sizeout(3) = sizeout_orig(3) + 2*(T(3));

HD = NaN(sizeout(1),sizeout(2),sizeout(3));
HD_temp = HD; HD_temp(1:size(hd,1),1:size(hd,2),1:size(hd,3)) = hd; hd = HD_temp;


Grid_Sim = grid;

TI2 = TI(1:end-T(1),1:end-T(2),1:end-T(3));


temp = ones(T(1),T(2), OL(3));
CCb = convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % top

temp = ones(OL(1), T(2), T(3));
CCf = convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % front

temp = ones(T(1),OL(2),T(3));
CCl = convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % left

temp = ones(OL(1),T(2)-OL(2), T(3)-OL(3));
CCfss = convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % front_small (k,j) for i,j,k

temp = ones(OL(1),T(2), T(3)-OL(3));
CCfs = convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % front_small (k) for i,k

temp=ones(T(1), OL(2),T(3)-OL(3));
CCls=convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % left_small (k)

temp=ones(T(1)-OL(1), OL(2),T(3));
CClm=convnfft(TI2.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % left_medium (i) for i,j


cntr = 0;
a = numel(1:T(1)-OL(1):sizeout(1)-T(1)+1);
b = numel(1:T(2)-OL(2):sizeout(2)-T(2)+1);
c = numel(1:T(3)-OL(3):sizeout(3)-T(3)+1);
LOC = NaN(a*b*c,2);


          
          cntr = cntr+1;
          
          selectrows = i:min(size(hd,1),i+CT(1)*T(1)-1);
          selectcols = j:min(size(hd,2),j+CT(2)*T(2)-1);
          selecthigh = k:min(size(hd,3),k+CT(3)*T(3)-1);
      
          hd_dev=hd(i:i+T(1)-1,j:j+T(2)-1,k:k+T(3)-1);
          hd_dev=hd_dev(:);
          hd_indicator=sum(isfinite(hd_dev));
          hd_index=find(~isnan(hd_dev));
          
          HD_dev=hd(selectrows, selectcols, selecthigh);
          HD_dev(1:T(1),1:T(2),1:T(3)) = NaN;
          HD_indicator=sum(isfinite(HD_dev(:)));
          
          dev_init = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+T(3)-1);
          
          
          if (i > 1) && (j > 1) && (k > 1),
              % top
              shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
              CC = CCb - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);

              % left
              shared = Grid_Sim(i:i+T(1)-1,j:j+OL(2)-1,k+OL(3):k+T(3)-1);
              CC2 = CCls - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC2(T(1):end-T(1)+1, OL(2):end-T(2)+1,T(3):end-T(3)+OL(3)+1);
        
              % front
              shared=Grid_Sim(i:i+OL(1)-1, j+OL(2):j+T(2)-1, k+OL(3):k+T(3)-1);
              CC3=CCfss-2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC3(OL(1):end-T(1)+1, T(2):end-T(2)+OL(2)+1, T(3):end-T(3)+OL(3)+1);

              % disp([num2str([i j k]),'---i j k---', 'ijk'])
          
          elseif i > 1 && j > 1
              % left
              shared = Grid_Sim(i+OL(1):i+T(1)-1,j:j+OL(2)-1,k:k+T(3)-1);
              CC = CClm - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+OL(1)+1, OL(2):end-T(2)+1,T(3):end-T(3)+1);
        
              %  front
              shared = Grid_Sim(i:i+OL(1)-1, j:j+T(2)-1, k:k+T(3)-1);
              CC2 = CCf-2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC2(OL(1):end-T(1)+1, T(2):end-T(2)+1, T(3):end-T(3)+1);
              
              % disp([num2str([i j k]),'---i j k---', 'ij'])

        elseif j > 1 && k > 1
           shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
           CC = CCb - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
           CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);

           % left
           shared = Grid_Sim(i:i+T(1)-1,j:j+OL(2)-1,k+OL(3):k+T(3)-1);
           CC2 = CCls - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
           CC = CC + CC2(T(1):end-T(1)+1, OL(2):end-T(2)+1,T(3):end-T(3)+OL(3)+1);
        
           % disp([num2str([i j k]),'---i j k---', 'jk'])
       elseif i > 1 && k > 1
           shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
           CC = CCb - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
           CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);
        
           shared = Grid_Sim(i:i+OL(1)-1, j:j+T(2)-1, k+OL(3):k+T(3)-1);
           CC2 = CCfs-2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
           CC = CC + CC2(OL(1):end-T(1)+1, T(2):end-T(2)+1, T(3):end-T(3)+OL(3)+1);
           % disp([num2str([i j k]),'---i j k---', 'ik'])
    
       elseif i > 1 
           shared=Grid_Sim(i:i+OL(1)-1, j:j+T(2)-1, k:k+T(3)-1);
           CC=CCf-2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
           CC = CC(OL(1):end-T(1)+1, T(2):end-T(2)+1, T(3):end-T(3)+1);
           % disp([num2str([i j k]),'---i j k---', 'i'])
           
       elseif j > 1
           shared = Grid_Sim(i:i+T(1)-1,j:j+OL(2)-1,k:k+T(3)-1);
           CC = CCl - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
           CC = CC(T(1):end-T(1)+1,OL(2):end-T(2)+1,T(3):end-T(3)+1);
           % disp([num2str([i j k]),'---i j k---', 'j'])
           
       elseif k > 1
           shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
           CC = CCb - 2 * convnfft(TI2,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
           CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);
           % disp([num2str([i j k]),'---i j k---', 'k'])
           
          else
              CC = rand(size(TI2,1)-T(1),size(TI2,2)-T(2),size(TI2,3)-T(3));
              % disp([num2str([i j k]),'---i j k---', 'else'])
          end;
          
          
          if hd_indicator==0 && HD_indicator==0
              [~,loc] = sort(CC(:));
              [x,y,z] = ind2sub(size(CC),loc(1:cand,1));
              if fc~=0
                  [c, ~] = hist_cat(TI, Grid_Sim, T, OL, fc, ibest, jbest, i, j);
              else
                  c = ceil(rand * length(x));
              end;
              pos = [x(c) y(c) z(c)];
              X_final = pos(1); Y_final = pos(2); Z_final = pos(3);
          else
              [~, loc] = sort(CC(:));
%               loc = loc(1:cand,1);
              [loc(:,2),loc(:,3),loc(:,4)] = ind2sub(size(CC),loc(:,1));
              now=0; difh=1E+5; difH=1E+5; DifT=1E+5;
              while (DifT~=0) && (now+1<=ceil(prop*size(loc,1))),
                  now = now+1;
                  x = loc(now,2);
                  y = loc(now,3);
                  z = loc(now,4);
                  target = TI2(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                  X = x:min(size(TI,1),x+CT(1)*T(1)-1);
                  Y = y:min(size(TI,2),y+CT(2)*T(2)-1);
                  Z = z:min(size(TI,3),z+CT(3)*T(3)-1);
                  TARGET = TI(X,Y,Z);
                  HD_dev2 = HD_dev(1:min(size(TARGET,1),size(HD_dev,1)),1:min(size(TARGET,2),size(HD_dev,2)),...
                      1:min(size(TARGET,3),size(HD_dev,3)));
                  HD_dev2(1:T(1),1:T(2),1:T(3)) = NaN;
                  HD_dev2 = reshape(HD_dev2,1,numel(HD_dev2));
                  HD_index2 = find(~isnan(HD_dev2));
                  Difh = sum(abs(hd_dev(hd_index)-target(hd_index)));
                  DifH = sum(abs(HD_dev2(HD_index2)-TARGET(HD_index2)));
                  DifT = sum(Difh + DifH);
                  
                  if Difh <= difh
                      difh = Difh;
                      X_final = x; Y_final = y; Z_final = z;
                      if Difh < difh
                          difH = 1E+5;
                      end;
                      if DifH < difH
                          difH = DifH;
                          X_final = x; Y_final = y; Z_final = z;
                      end;
                  end;
              end;
          end;
          
          LOC (cntr,1)= X_final; LOC (cntr,2)= Y_final; LOC (cntr,3)= Z_final;

          Target = TI2(X_final:X_final+T(1)-1,Y_final:Y_final+T(2)-1,Z_final:Z_final+T(3)-1);
          
          
          M = mincut3D_func(Target, Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+T(3)-1), T, OL, i, j, k );
          
          dev_index = isnan(dev_init);
          dev_init (dev_index) = 0;
          
          Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+T(3)-1) = combine(dev_init,Target,M);
      
