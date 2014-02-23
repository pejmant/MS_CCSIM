
% -------------------------------------------------------------------------
% Reference: Tahmasebi, P., Sahimi, M., Caers, J., 2013. 
% MS-CCSIM: accelerating pattern-based geostatistical simulation of 
% categorical variables using a multi-scale search in Fourier space,
% Computers & Geosciences, 


% Author: Pejman Tahmasebi
% E-mail: pejman@stanford.edu
% Stanford Center for reservoir Forecasting, Energy Resources Engineering 
% Department,Stanford University.
% -----------------------------------------------------------------------*/


function [Grid_Sim, LOC2] = CCSIM_3D_MS1(TI, hd, LOC1, T, OL, fc, rad)

%% Input Parameters
% - TI: Training image
% - hd: hard data matrix
% - LOC1: Location of the previous coarser patterns
% - T: Template size
% - OL: Overlap size
% - rad: Radius of pattern searching

%% Output Parameters
% - Grid_Sim: Simulation grid for output

%% ---------------------------------------------------------------------------------------- 

sizeout = size(hd);

Grid_Sim = zeros(sizeout);

cntr = 0;
a = numel(1:T(1)-OL(1):sizeout(1)-T(1)+1);
b = numel(1:T(2)-OL(2):sizeout(2)-T(2)+1);
c = numel(1:T(3)-OL(3):sizeout(3)-T(3)+1);
LOC2 = NaN(a*b*c,2);
index1 = reshape(1:size(TI,1)*size(TI,2)*size(TI,3),size(TI,1),size(TI,2),size(TI,3));


for i=[1:T(1)-OL(1):sizeout(1)-T(1), sizeout(1)-T(1)+1],
  for j=[1:T(2)-OL(2):sizeout(2)-T(2), sizeout(2)-T(2)+1],
      for k=[1:T(3)-OL(3):sizeout(3)-T(3), sizeout(3)-T(3)+1],
          
          cntr = cntr+1;
          hd_dev = hd(i:i+T(1)-1,j:j+T(2)-1,k:k+T(3)-1);
          hd_dev = hd_dev(:);
          hd_indicator = sum(isfinite(hd_dev));
          hd_index = find(~isnan(hd_dev));
     
     
          LOC1_x = 2*LOC1(cntr,1)-1; LOC1_y = 2*LOC1(cntr,2)-1; LOC1_z = LOC1(cntr,3);
          selectrows = max(1,LOC1_x-rad*T(1)):min(size(TI,1),LOC1_x+rad*T(1)+2);
          selectcols = max(1,LOC1_y-rad*T(2)):min(size(TI,2),LOC1_y+rad*T(2)+2);
          selecthight = max(1,LOC1_z-rad*T(3)):min(size(TI,3),LOC1_z+rad*T(3));
          domain = TI(selectrows, selectcols, selecthight);
          index_section = index1(selectrows, selectcols, selecthight);
          
          
          if (i > 1) && (j > 1) && (k > 1),
              temp = ones(T(1),T(2), OL(3));
              CCb = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % bot

              temp = ones(T(1), OL(2),T(3)-OL(3));
              CCls = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % side_left_small

              temp = ones(OL(1),T(2)-OL(2), T(3)-OL(3));
              CCfss = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % side_front_small
       
       
              shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
              CC = CCb - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);
              
              shared = Grid_Sim(i:i+T(1)-1,j:j+OL(2)-1,k+OL(3):k+T(3)-1);
              CC2 = CCls - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC2(T(1):end-T(1)+1, OL(2):end-T(2)+1,T(3):end-T(3)+OL(3)+1);
       
              shared=Grid_Sim(i:i+OL(1)-1, j+OL(2):j+T(2)-1, k+OL(3):k+T(3)-1);
              CC3=CCfss-2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC3(OL(1):end-T(1)+1, T(2):end-T(2)+OL(2)+1, T(3):end-T(3)+OL(3)+1);
              
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:));
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end
                  
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          
          elseif i > 1 && j > 1
              temp = ones(OL(1), T(2), T(3));
              CCf = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % front

              temp=ones(T(1)-OL(1), OL(2),T(3));
              CClm=convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % left_medium (i) for i,j

              shared = Grid_Sim(i+OL(1):i+T(1)-1,j:j+OL(2)-1,k:k+T(3)-1);
              CC = CClm - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+OL(1)+1, OL(2):end-T(2)+1,T(3):end-T(3)+1);
       
              shared=Grid_Sim(i:i+OL(1)-1, j:j+T(2)-1, k:k+T(3)-1);
              CC2=CCf-2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC2(OL(1):end-T(1)+1, T(2):end-T(2)+1, T(3):end-T(3)+1);
     
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:));
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          
          elseif j > 1 && k > 1
              temp = ones(T(1),T(2), OL(3));
              CCb = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % top

              temp=ones(T(1), OL(2),T(3)-OL(3));
              CCls=convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % left_small (k)       
       
              shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
              CC = CCb - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);

              shared = Grid_Sim(i:i+T(1)-1,j:j+OL(2)-1,k+OL(3):k+T(3)-1);
              CC2 = CCls - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC + CC2(T(1):end-T(1)+1, OL(2):end-T(2)+1,T(3):end-T(3)+OL(3)+1);
              
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);           
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:));
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          
          elseif i > 1 && k > 1
              temp = ones(T(1),T(2), OL(3));
              CCb = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % top

              temp = ones(OL(1),T(2), T(3)-OL(3));
              CCfs = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % front_small (k) for i,k
       
              shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
              CC = CCb - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);
       
              shared=Grid_Sim(i:i+OL(1)-1, j:j+T(2)-1, k+OL(3):k+T(3)-1);
              CC2=CCfs-2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC= CC + CC2(OL(1):end-T(1)+1, T(2):end-T(2)+1, T(3):end-T(3)+OL(3)+1);
              
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);           
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:));
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          
          elseif i > 1
              temp = ones(OL(1), T(2), T(3));
              CCf = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % front
       
              shared=Grid_Sim(i:i+OL(1)-1, j:j+T(2)-1, k:k+T(3)-1);
              CC=CCf-2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
              CC = CC(OL(1):end-T(1)+1, T(2):end-T(2)+1, T(3):end-T(3)+1);
              
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);           
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:)); 
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          
          elseif j > 1
              temp = ones(T(1),OL(2),T(3));
              CCl = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % left
       
              shared = Grid_Sim(i:i+T(1)-1,j:j+OL(2)-1,k:k+T(3)-1);
              CC = CCl - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1)) + sum(shared(:).^2);
     
              CC = CC(T(1):end-T(1)+1,OL(2):end-T(2)+1,T(3):end-T(3)+1);
              
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);           
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:));
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          elseif k > 1
              temp = ones(T(1),T(2), OL(3));
              CCb = convnfft(domain.^2, temp(end:-1:1,end:-1:1,end:-1:1)); % top
       
              shared = Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+OL(3)-1);
              CC = CCb - 2 * convnfft(domain,shared(end:-1:1,end:-1:1,end:-1:1))+ sum(shared(:).^2);
              CC = CC(T(1):end-T(1)+1,T(2):end-T(2)+1,OL(3):end-T(3)+1);
              
              if hd_indicator==0
                  [x,y,z] = ind2sub(size(CC),find(CC==min(CC(:))));
                  if fc~=0
                      [c, ~] = hist_3D_cat(TI, Grid_Sim, T, OL, fc, x, y, z, i, j, k);
                  else
                      c = ceil(rand * length(x));
                  end;
                  target_final = domain(x(c):x(c)+T(1)-1, y(c):y(c)+T(2)-1, z(c):z(c)+T(3)-1);
                  index_intrst = index_section(x(c),y(c),z(c));
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);           
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              else %find a pattern in patdb base on its HD and CC
                  [~, loc]=sort(CC(:));
                  % Identify the location of scores
                  [loc(:,2),loc(:,3),loc(:,4)]=ind2sub(size(CC),loc(:));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(now,2);
                      y=loc(now,3);
                      z=loc(now,4);
                      target=domain(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          
          
          else
              if hd_indicator==0
                  target_final = TI(LOC1_x:LOC1_x+T(1)-1,LOC1_y:LOC1_y+T(2)-1,LOC1_z:LOC1_z+T(3)-1);
                  LOC2 (cntr,1)=LOC1_x; LOC2 (cntr,2)=LOC1_y; LOC2 (cntr,3)=LOC1_z;
              else
                  v = NaN(size(domain));
                  v(end-T(1):end,:,:) = 1; v(:,end-T(2):end,:) = 1; v(:,:,end-T(3):end) = 1;
                  [v_x,v_y,v_z] = ind2sub(size(v),find(isnan(v)==1));
                  loc = v_x; loc(:,2)=v_y; loc(:,3)=v_z;
                  R = randperm(numel(v_x));
                  now=0; dif=1E+5;
                  while (dif~=0 && now+1<=size(loc,1)),
                      now=now+1;
                      x=loc(R(now),1);
                      y=loc(R(now),2);
                      z=loc(R(now),3);
                      target=TI(x:x+T(1)-1,y:y+T(2)-1,z:z+T(3)-1);
                      D = sum(abs(hd_dev(hd_index)-target(hd_index)));
                      if D < dif
                          dif = D;
                          target_final = target;
                      end;
                  end;
                  index_intrst = index_section(x,y,z);
                  [xx,yy,zz]=ind2sub(size(TI),index_intrst);             
                  LOC2 (cntr,1)=xx; LOC2 (cntr,2)=yy; LOC2 (cntr,3)=zz;
              end;
          end;
          
          M = mincut3D_func( target_final, Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+T(3)-1), T, OL, i, j, k );
          Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+T(3)-1) = combine(Grid_Sim(i:i+T(1)-1, j:j+T(2)-1, k:k+T(3)-1),target_final, M);
      end;
  end;
end;





