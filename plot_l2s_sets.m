%{
[sergio@c24-52 sarta_utilities]$ ls -lt /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/*/AIRS*_1.mat
-rw-rw-r-- 1 sergio pi_sergio 118679766 Aug  5 06:53 /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/FWOsun/AIRS_L1C_R49_allPaths_14a_1.mat
-rw-rw-r-- 1 sergio pi_sergio 148044147 Aug  5 06:53 /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/FCOW/AIRS_L1C_R49_allPaths_14a_1.mat
-rw-rw-r-- 1 sergio pi_sergio  43321506 Aug  5 06:53 /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/wvFMW/AIRS_L1C_R49_allPaths_14a_1.mat
-rw-rw-r-- 1 sergio pi_sergio  58265087 Aug  5 06:53 /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/FOW/AIRS_L1C_R49_allPaths_14a_1.mat
-rw-rw-r-- 1 sergio pi_sergio  58265087 Aug  5 06:53 /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/FWO/AIRS_L1C_R49_allPaths_14a_1.mat

FWOsun : {'F/convolved_kcarta_F_'}    {'FO/convolved_kcarta_FO_'}    {'FWO/convolved_kcarta_FWO_'}    {'FWOP/convolved_kcarta_FWOP_'}
FCOW   : {'cobandF/convolved_k...'}    {'cobandFC/convolved_...'}    {'cobandFCO/convolved...'}    {'cobandFCOW/convolve...'}    {'FWOP/convolved_kcar...'}
wvFMW  : {'wvbandF/convolved_kcarta_wvbandF_'}    {'wvbandFM/convolved_kcarta_FO_'}    {'wvbandFMW/convolved_kcarta_FWO_'}
FOW    : {'F/convolved_kcarta_F_'}    {'FO/convolved_kcarta_FO_'}    {'FWO/convolved_kcarta_FWO_'}    {'FWOP/convolved_kcarta_FWOP_'}
FWO    : {'F/convolved_kcarta_F_'}    {'FO/convolved_kcarta_FO_'}    {'FWO/convolved_kcarta_FWO_'}    {'FWOP/convolved_kcarta_FWOP_'}

Then use the very simple reader for the .mat files :
  plot_l2s_sets.m

See ftc_dev_how_to_create_sarta_ftc.pdf,how_to_create_sarta_ftc.pdf
SARTA has 7 bands/sets
  Sets 1,2     : FWO,FOW      : LW band : needs F,FO,FWO,FWOP
  Set 3        : wvFMW        : MW band : need F,FM,FMW
  Sets 4,5,6,7 : CFOW, FWOsun : SW band : need F,FC,FCO,FCOW,FCOWP

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% script to look at l2s tranmittances made by doall_wrtconvdat_generic
%% script to look at l2s tranmittances made by doall_wrtconvdat_generic
%% script to look at l2s tranmittances made by doall_wrtconvdat_generic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%load /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/wvFMW/AIRS_L1C_R49_allPaths_14a_1.mat
load /home/sergio/nogit/ftcprod/prod_2026/airs_l1c/aug2026_regr49_aircraft_12km_EXPT/FWO/AIRS_L1C_R49_allPaths_14a_1.mat

mfiles

fprintf(1,'number of angs = %2i number of mfiles = %2i \n',length(secangs),length(mfiles))
disp(' ')
printarray(size(L2Sdata.ctrans),'size(L2Sdata.ctrans)')

%%%%%%%%%%%%%%%%%%%%%%%%%
disp('plotting first mfile set v1')
for iang = 1 : length(secangs)
  ibreak = 1;
  offsetang   = iang-1;
  offsetbreak = ibreak-1;
  
  indoffset = offsetbreak * length(secangs) + offsetang;
  indx = (1:100) + indoffset*100;  
  imagesc(L2Sdata.f,1:100,L2Sdata.ctrans(indx,:)); colorbar; colormap jet
  title(['v1 angsles ' num2str(iang) ' of ' num2str(length(secangs))])
  disp('ret to continue'); pause  
  %pause(1)
end

disp('plotting first mfile set v2')
for ibreak = 1 : length(mfiles)
  iang = 1;
  offsetang   = iang-1;
  offsetbreak = ibreak-1;
  
  indoffset = offsetbreak * length(secangs) + offsetang;
  indx = (1:100) + indoffset*100;  
  imagesc(L2Sdata.f,1:100,L2Sdata.ctrans(indx,:)); colorbar; colormap jet
  title(['v2 ' num2str(ibreak)])  
  title(num2str(ibreak))
  title(['v2 breakouts ' num2str(iang) ' of ' num2str(length(mfiles)) ' : ' mfiles{ibreak}])  
  disp('ret to continue'); pause    
  %pause(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%

iang   = input('Enter angle to plot : (-1 to end) ');
ibreak = input('Enter mfile to plot : (-1 to end) ');

figure(1); clf
while iang > 0 & iang <= length(secangs) & ibreak > 0 & ibreak <= length(mfiles) 
  ind = 1;
  offsetang   = iang-1;
  offsetbreak = ibreak-1;
  
  indoffset = offsetbreak * length(secangs) + offsetang;
  %% if iang = 1,ibreak = 1 this works out to 0;
  %% if iang = 2,ibreak = 1 this works out to 1;
  %% if iang = 1,ibreak = 2 this works out to length(scanangs) so goes to next breakout set;    
  
  indx = (1:100) + (ind-1)*100;
  indx = (1:100) + indoffset*100;  
  imagesc(L2Sdata.f,1:100,L2Sdata.ctrans(indx,:)); colorbar; colormap jet; title(mfiles{ibreak})

  disp(' ')
  iang   = input('Enter angle to plot : (-1 to end) ');
  ibreak = input('Enter mfile to plot : (-1 to end) ');
end
