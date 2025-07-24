function [iok] = write_trace_gas_channel_lists(csens)

%
% NB for CrIS the sorting of guard chanels comes after creating the
% fast coefficients using reorder_gurad_chans.m
%

csens=upper(csens)
if(~ismember(csens,{'AIRS_L1C','CRIS_LR','CRIS_HR','CHIRP','IASI'}))
  error('Invalid sensor choice')
  return;
end

switch csens
  case 'CRIS_LR'
    x=load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    dout = '/home/chepplew/data/sarta/tmp_chanlists/';
    crfn = 'cris_lr_g4_list_';
        
  case 'CRIS_HR'
  
  case 'CHIRP'
  
  case 'AIRS_L1C'
  
  case 'IASI'
  
end

%
% --------------------------------------------
% SO2 
so2.wnbnd = [1059 1245; 1301 1404; 2450 2530];   % cm-1. 3 bands ([1,2,3],:)
    
ix1 = find(x.vchan >= so2.wnbnd(1,1) & x.vchan <= so2.wnbnd(1,2));
ix2 = find(x.vchan >= so2.wnbnd(2,1) & x.vchan <= so2.wnbnd(2,2));
ix3 = find(x.vchan >= so2.wnbnd(3,1) & x.vchan <= so2.wnbnd(3,2));
ix = [ix1 ix2 ix3];
[~,ichan] = intersect(x.vchan, x.vchan(ix));

fnout = [dout crfn 'so2'];
FH = fopen(fnout,'w')
for i = 1:length(ix)
%  fprintf(FH,'   %3d   %8.3f\n', x.idchan(ix(i)), x.vchan(ix(i)) );
  fprintf(FH,'   %3d   %8.3f\n', ichan(i), x.vchan(ix(i)) );
end

% ------------------------------------------------
% HNO3
hno3.wnbnd = [740.00 775.00; 850.00 925.00; 1212.00 1356; 1670.00 1745.00 ];
ix1 = find(x.vchan >= hno3.wnbnd(1,1) & x.vchan <= hno3.wnbnd(1,2));
ix2 = find(x.vchan >= hno3.wnbnd(2,1) & x.vchan <= hno3.wnbnd(2,2));
ix3 = find(x.vchan >= hno3.wnbnd(3,1) & x.vchan <= hno3.wnbnd(3,2));
ix4 = find(x.vchan >= hno3.wnbnd(4,1) & x.vchan <= hno3.wnbnd(4,2));
ix = [ix1 ix2 ix3 ix4];
[~,ichan] = intersect(x.vchan, x.vchan(ix));

fnout = [dout crfn 'hno3'];
FH = fopen(fnout,'w')
for i = 1:length(ix)
%  fprintf(FH,'   %3d   %8.3f\n', x.idchan(ix(i)), x.vchan(ix(i)) );
  fprintf(FH,'   %3d   %8.3f\n', ichan(i), x.vchan(ix(i)) );
end

% ----------------------------------------------------
% N2O (247 chans)
n2o.wnbnd = [1126 1205; 1230 1328; 1838 1910; 2115 2560];
ix1 = find(x.vchan >= n2o.wnbnd(1,1) & x.vchan <= n2o.wnbnd(1,2));
ix2 = find(x.vchan >= n2o.wnbnd(2,1) & x.vchan <= n2o.wnbnd(2,2));
ix3 = find(x.vchan >= n2o.wnbnd(3,1) & x.vchan <= n2o.wnbnd(3,2));
ix4 = find(x.vchan >= n2o.wnbnd(4,1) & x.vchan <= n2o.wnbnd(4,2));
ix = [ix1 ix2 ix3 ix4];
[~,ichan] = intersect(x.vchan, x.vchan(ix));

fnout = [dout crfn 'n2o'];
FH = fopen(fnout,'w')
for i = 1:length(ix)
%  fprintf(FH,'   %3d   %8.3f\n', x.idchan(ix(i)), x.vchan(ix(i)) );
  fprintf(FH,'   %3d   %8.3f\n', ichan(i), x.vchan(ix(i)) );
end

% ----------------------------------------------------
% NH3 ( chans)
nh3.wnbnd = [750 1220; 1450 1800];
ix1 = find(x.vchan >= nh3.wnbnd(1,1) & x.vchan <= nh3.wnbnd(1,2));
ix2 = find(x.vchan >= nh3.wnbnd(2,1) & x.vchan <= nh3.wnbnd(2,2));
ix = [ix1 ix2];
[~,ichan] = intersect(x.vchan, x.vchan(ix));

fnout = [dout crfn 'nh3'];
FH = fopen(fnout,'w')
for i = 1:length(ix)
%  fprintf(FH,'   %3d   %8.3f\n', x.idchan(ix(i)), x.vchan(ix(i)) );
  fprintf(FH,'   %3d   %8.3f\n', ichan(i), x.vchan(ix(i)) );
end



