
% -------------------------------------------------------------------------
% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space
% functions, Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Stanford University.
% -----------------------------------------------------------------------*/


function [Grid_Sim, LOC2] = CCSIM_2D_MS1(TI, hd, LOC1, T, OL, rad)

sizeout = size(hd);
Grid_Sim = zeros(sizeout);

cntr = 0;
a = numel(1:T-OL:sizeout(1)-T+1);
b = numel(1:T-OL:sizeout(2)-T+1);
LOC2 = NaN(a*b,2);
index1 = reshape(1:size(TI,1)*size(TI,2),size(TI,1),size(TI,2));

% tic;
for i=[1:T-OL:sizeout(1)-T, sizeout(1)-T+1],
  for j=[1:T-OL:sizeout(2)-T, sizeout(2)-T+1],
      
      cntr = cntr+1;
      
      hd_dev = hd(i:i+T-1,j:j+T-1);
      hd_dev = reshape(hd_dev,1,size(hd_dev,1)*size(hd_dev,2));
      hd_indicator = sum(isfinite(hd_dev));
      hd_index = find(~isnan(hd_dev));
      
      LOC1_x = 2*LOC1(cntr,1)-1; LOC1_y = 2*LOC1(cntr,2)-1;
      siz = size(TI);
      selectrows = max(1,LOC1_x-rad*T):min(siz(1),LOC1_x+rad*T);
      selectcols = max(1,LOC1_y-rad*T):min(siz(2),LOC1_y+rad*T);
      domain = TI(selectrows, selectcols);
      index_section = index1(selectrows, selectcols);      
            

      if (i > 1) && (j > 1),
          
          temp = ones([OL T]);
          CCtop = xcorr2(domain.^2, temp);
          temp = ones([T-OL OL]);
          CCsidesmall = xcorr2(domain.^2, temp);
          
          shared = Grid_Sim(i:i+OL-1,j:j+T-1);
          CC = CCtop - 2 * xcorr2(domain, shared) + sum(shared(:).^2);
          CC = CC(OL:end-T+1,T:end-T+1);
          shared = Grid_Sim(i+OL:i+T-1,j:j+OL-1);
          CC2 = CCsidesmall - 2 * xcorr2(domain, shared) + sum(shared(:).^2);
          CC = CC + CC2(T:end-T+OL+1, OL:end-T+1);
               
      
      elseif i > 1
          
          temp = ones([OL T]);
          CCtop = xcorr2(domain.^2, temp);      
          
          shared = Grid_Sim(i:i+OL-1,j:j+T-1);
          CC = CCtop - 2 * xcorr2(domain, shared) + sum(shared(:).^2);

          CC = CC(OL:end-T+1,T:end-T+1);

    elseif j > 1
          temp = ones([T OL]);
          CCside = xcorr2(domain.^2, temp);
        
          shared = Grid_Sim(i:i+T-1,j:j+OL-1);
          CC = CCside - 2 * xcorr2(domain, shared) + sum(shared(:).^2);
      
          CC = CC(T:end-T+1,OL:end-T+1);        
          
          
      else
          CC = zeros(size(domain));
          CC(end-T+1:end,:) = 1; 
          CC(:,end-T+1:end) = 1;
      end;
      
      if hd_indicator==0 && (i||j)~=1
          [ibest, jbest] = find(CC == min(CC(:)));
          c = ceil(rand * length(ibest));
          pos = [ibest(c) jbest(c)];
          target_final = domain(pos(1):pos(1)+T-1,pos(2):pos(2)+T-1);
          index_intrst = index_section(ibest(c),jbest(c));
          [xx,yy]=ind2sub(size(TI),index_intrst);
          LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy;
      elseif hd_indicator==0 && (i&&j)==1
          target_final = TI(LOC1_x:LOC1_x+T-1,LOC1_y:LOC1_y+T-1);
          LOC2 (cntr,1)=LOC1_x; LOC2 (cntr,2)=LOC1_y;
      else %find a pattern in patdb base on its HD and CC
          [~, loc]=sort(CC(:));
          % Identify the location of scores
          [loc(:,2),loc(:,3)]=ind2sub(size(CC),loc(:));
          now=0; dif=1E+5;
          while (dif~=0 && now+1<=size(loc,1)),
              now=now+1;
              x=loc(now,2);
              y=loc(now,3);
              target=domain(x:x+T-1,y:y+T-1);
              D = sum(abs(hd_dev(hd_index)-target(hd_index)));
              if D < dif
                  dif = D;
                  target_final = target;
              end;
          end;
          index_intrst = index_section(x,y);
          [xx,yy]=ind2sub(size(TI),index_intrst);
          LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy;
      end;
      
      M = mincut_func( target_final, Grid_Sim(i:i+T-1, j:j+T-1), T, OL, i, j );
      Grid_Sim(i:i+T-1,j:j+T-1) = combine_2D(Grid_Sim(i:i+T-1, j:j+T-1),target_final, M);
  end;
end;














