
% -------------------------------------------------------------------------
% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences,  


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Stanford University.
% -----------------------------------------------------------------------*/


function [Grid_Sim, LOC] = CCSIM_2D_MS2(TI, hd, T, OL, CT, fc, prop, cand)

sizeout = size(hd);
Grid_Sim = NaN(sizeout);

TI2 = TI(1:end-T,1:end-T);

temp = ones([OL T]);
CCtop = xcorr2(TI2.^2, temp);
temp = ones([T OL]);
CCside = xcorr2(TI2.^2, temp);
temp = ones([T-OL OL]);
CCsidesmall = xcorr2(TI2.^2, temp);

cntr = 0;
a = numel(1:T-OL:Grid_Sim(1)-T+1);
b = numel(1:T-OL:Grid_Sim(2)-T+1);
LOC = NaN(a*b,2);


for i=[1:T-OL:sizeout(1)-T, sizeout(1)-T+1],
  for j=[1:T-OL:sizeout(2)-T, sizeout(2)-T+1],
      
      cntr = cntr+1;
      
      selectrows = i:min(size(hd,1),i+CT(1)*T-1);
      selectcols = j:min(size(hd,2),j+CT(2)*T-1);
          
      hd_dev = hd(i:i+T-1,j:j+T-1);
      hd_dev = reshape(hd_dev,1,size(hd_dev,1)*size(hd_dev,2));
      hd_indicator = sum(isfinite(hd_dev));
      hd_index = find(~isnan(hd_dev));
      
      HD_dev = hd(selectrows, selectcols);
      HD_dev(1:T,1:T) = NaN;
      HD_indicator = sum(isfinite(HD_dev(:)));
      
      dev_init = Grid_Sim(i:i+T-1, j:j+T-1);

      if (i > 1) && (j > 1),
          shared = Grid_Sim(i:i+OL-1,j:j+T-1);
          CC = CCtop - 2 * xcorr2(TI2, shared) + sum(shared(:).^2);
          CC = CC(OL:end-T+1,T:end-T+1);

          shared = Grid_Sim(i+OL:i+T-1,j:j+OL-1);
          CC2 = CCsidesmall - 2 * xcorr2(TI2, shared) + sum(shared(:).^2);
          CC = CC + CC2(T:end-T+OL+1, OL:end-T+1);
    
      elseif i > 1
          shared = Grid_Sim(i:i+OL-1,j:j+T-1);
          CC = CCtop - 2 * xcorr2(TI2, shared) + sum(shared(:).^2);
          CC = CC(OL:end-T+1,T:end-T+1);

      elseif j > 1
        shared = Grid_Sim(i:i+T-1,j:j+OL-1);
        CC = CCside - 2 * xcorr2(TI2, shared) + sum(shared(:).^2);
        CC = CC(T:end-T+1,OL:end-T+1);
      
      else
          CC = rand(size(TI2,1)-T,size(TI2,2)-T);
      end;
      
      if hd_indicator==0 && HD_indicator==0
          [~,loc] = sort(CC(:));
          [ibest, jbest] = ind2sub(size(CC),loc(1:cand,1));
          if fc~=0
              [c, ~] = hist_cat(TI, Grid_Sim, T, OL, fc, ibest, jbest, i, j);
          else
              c = ceil(rand * length(ibest));
          end;
          pos = [ibest(c) jbest(c)];
          X_final = pos(1); Y_final = pos(2);
%           LOC (cntr,1) = pos(1); LOC (cntr,2) = pos(2);
      else
          [~, loc] = sort(CC(:));
          [loc(:,2),loc(:,3)] = ind2sub(size(CC),loc(:));
          now=0; difh=1E+5; difH=1E+5; DifT=1E+5; 
              
          while (DifT~=0) && (now+1<=ceil(prop*size(loc,1))),
              now = now+1;
              x = loc(now,2);
              y = loc(now,3);
              target = TI2(x:x+T-1,y:y+T-1);
              X = x:min(size(TI,1),x+CT(1)*T-1);
              Y = y:min(size(TI,2),y+CT(2)*T-1);
              TARGET = TI(X, Y);
              HD_dev2 = HD_dev(1:min(size(TARGET,1),size(HD_dev,1)),1:min(size(TARGET,2),size(HD_dev,2)));
              HD_dev2(1:T,1:T) = NaN;
              HD_dev2 = reshape(HD_dev2,1,size(HD_dev2,1)*size(HD_dev2,2));
              HD_index2 = find(~isnan(HD_dev2));
              Difh = sum(abs(hd_dev(hd_index)-target(hd_index)));
              DifH = sum(abs(HD_dev2(HD_index2)-TARGET(HD_index2)));
              DifT = sum(Difh + DifH);
              
              if Difh <= difh
                  difh = Difh;
                  X_final = x; Y_final = y;
                  if Difh < difh
                      difH = 1E+5;
                  end;
                  if DifH < difH
                      difH = DifH;
                      X_final = x; Y_final = y;
                  end;
              end;
              
          end
      end;
      
      LOC (cntr,1)= X_final; LOC (cntr,2)= Y_final;
      
      Target = TI2(X_final:X_final+T-1,Y_final:Y_final+T-1);
      
      M = mincut_func(Target, Grid_Sim(i:i+T-1, j:j+T-1), T, OL, i, j);
      
      dev_index = isnan(dev_init);
      dev_init (dev_index) = 0;
      
      Grid_Sim(i:i+T-1,j:j+T-1) = combine_2D(dev_init,Target, M);
  end
end;