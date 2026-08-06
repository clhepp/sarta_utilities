% function make_refprof.m
%
% Workflow:
%  1. load the LAYER data for the regression set. The 49th is the us.std used
%     as the reference in SARTA.
%  2. Copy the header from an older refprof file, update if needed.
%  3. compute layer thickness and prof.plays values as needed.
%  4. Write the ascii-text reference profile file for SARTA
%  5. write the MAT reference profile file for water continuum modeling.
%
% Jan 2025   CLH : edit for airs_oco2_pbl       layering
% Aug 2026   SSM : edit for airs 12 km aircraft layering
%
% see ~/git/ftc_dev/Readme_ftc_notes_dryrun49Profiles_AIRS_aircraft_12km for example usage/need

addpath /home/sergio/git/matlabcode/matlibSergio/matlib2025/h4tools

iVers = 1;  %% CLH PBL file
iVers = 2;  %% SSM 12 km aircraft file
iVers = 3;  %% SSM 20 km aircraft file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iVers == 1
  %  Output path and file
  
  sav.dir    = '/home/chepplew/data/sarta/prod_2025/generic/';
  sav.fn_txt = 'refprof_400ppm_pbl.txt';
  sav.fn_mat = 'refprof_400ppm_pbl.mat';
  
  % ----------------------------------------------------------------------
  %    Update the reference profile
  % ----------------------------------------------------------------------
  % get the 49 regression input profiles 
  if( strcmp(sens,'CRIS') )
    fnr49 = '/home/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm_2235.op.rtp';
  end
  if( strcmp(sens,'AIRS') )
    fnr49 = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr49_1100_400ppm.op.rtp';
  end
  
  % load airs_oco2_pbl LAYERS file
  fnr49 = ['/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/RUN_KCARTA/' ... 
           'REGR49_400ppm_H2020_Jan2025_PBL_AIRS2834_3CrIS_IASI/regr49_pbl.op.rtp'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%

elseif iVers == 2
  %% input profile
  fnr49 = '/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1100_with_co2_400ppm_9gases_unitemiss_aircraft_12km.op.rtp';

  %  Output path and file
  sav.dir    = '/home/sergio/git/ftc_dev/Refs/';
  sav.fn_txt = 'refprof_400ppm_aircraft_12km.txt';
  sav.fn_mat = 'refprof_400ppm_aircraft_12km.mat';

elseif iVers == 3
  %% input profile
  fnr49 = '/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1100_with_co2_400ppm_9gases_unitemiss_aircraft_20km.op.rtp';

  %  Output path and file
  sav.dir    = '/home/sergio/git/ftc_dev/Refs/';
  sav.fn_txt = 'refprof_400ppm_aircraft_20km.txt';
  sav.fn_mat = 'refprof_400ppm_aircraft_20km.mat';

else
  error('unknown iVers')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[head hattr prof pattr] = rtpread(fnr49);

% Use the reference file used by Scott for header lines
%reffn  = '/asl/data/sarta_database/Data_IASI_sep08/Coef/profref_trace385';
%reffn  = '/home/chepplew/projects/sarta/cris_hr/profref_trace385';  % extra comment line
%reffn = '/home/chepplew/data/sarta/prod_2019/iasi/dec2018/dbase/Coef/refprof_nh3';
reffn = '/home/sergio/git/ftc_dev/Refs/refprof_nh3';

inFH = fopen(reffn,'r');
A    = textscan(inFH,'%d %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','!');
fclose(inFH);
A = importdata(reffn);

% calculated layer thickness from regr49
for i=1:100 
 thick(i) = prof.palts(i,49) - prof.palts(i+1,49); 
end

% compare Sergio's 49th w/ Scott's:
AVOG    = 6.02214199E+26;
STDATM  = 1013.25;
figure(1); clf; plot(A.data(:,4),A.data(:,1),'o-',prof.plevs(100:-1:1,49)/1013.25,[1:100],'+-'); grid on; xlabel('P [atm]');       ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 100])
figure(2); clf; plot(A.data(:,3),A.data(:,1),'o-',thick(100:-1:1),[1:100],'+-');                 grid on; xlabel('lay thick [m]'); ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 100])
figure(3); clf; plot(A.data(:,5),A.data(:,1),'o-',prof.ptemp(101:-1:1,49),[1:101],'+-');         grid on; xlabel('T [K]');         ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 101])
figure(4); clf; plot(A.data(:,2),A.data(:,1),'o-',prof.palts(101:-1:1,49),[1:101],'+-');         grid on; xlabel('lay alt [m]');   ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 101])
figure(5); clf; plot(A.data(:,6),A.data(:,1),'o-',prof.gas_2(100:-1:1,49)/AVOG,[1:100],'+-');    grid on; xlabel('CO2 molec/cm2'); ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 100])

figure(1); clf; plot(A.data(:,4),A.data(:,1),'o-',prof.plevs(100:-1:1,49)/1013.25,[1:100],'+-'); grid on; xlabel('P [atm]');       ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 100])
figure(2); clf; plot(A.data(:,3),A.data(:,1),'o-',thick(100:-1:1),[1:100],'+-');                 grid on; xlabel('lay thick [m]'); ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 100])
figure(3); clf; plot(A.data(:,5),A.data(:,1),'o-',prof.ptemp(101:-1:1,49),[1:101],'+-');         grid on; xlabel('T [K]');         ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 101])
figure(4); clf; plot(A.data(:,2),A.data(:,1),'o-',prof.palts(101:-1:1,49),[1:101],'+-');         grid on; xlabel('lay alt [m]');   ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 101])
figure(5); clf; plot(A.data(:,6),A.data(:,1),'o-',prof.gas_2(100:-1:1,49)/AVOG,[1:100],'+-');    grid on; xlabel('CO2 molec/cm2'); ylabel('layer'); legend('US STD','new prof','location','best'); ylim([1 100])

% -----------------------------------------------------------------------
% convert pressure on levels to pressure of layer between the two levels.

pN = prof.plevs(1:100,:)-prof.plevs(2:101,:);
pD = log(prof.plevs(1:100,:) ./ prof.plevs(2:101,:));
prof.plays = zeros(size(prof.plevs));
prof.plays(1:100,:) = pN ./ pD;

%%! no altitude(m)  thickness(m) pressure(atm) temp(K) fixed=CO2     H2O        O3         CO         CH4        SO2       HNO3         N2O       NH3
%%   1     2           3           4             5       6            7         8           9         10          11        12           13       14
figure(6);clf;semilogy(A.data(:,6),A.data(:,4)*STDATM,'o-',prof.gas_2(1:100,49)/AVOG,prof.plays(1:100,49),'+-');       grid on; xlabel('CO2 amt'); ylabel('plays [mb]'); legend('US STD','new prof','location','best')
  grid on;set(gca,'YDir','Reverse'); ylim([5e-3 1100])
%figure(7);clf;semilogy(A.data(:,6),prof.plays(100:-1:1),'o-',prof.gas_2(1:100,49)/AVOG,prof.plays(1:100,49),'+-'); grid on; xlabel('CO2 amt'); ylabel('plays [mb]'); legend('US STD','new prof','location','best')
%  grid on;set(gca,'YDir','Reverse'); ylim([5e-3 1100])

figure(7);clf;semilogy(A.data(:,7),A.data(:,4)*STDATM,'o-',prof.gas_1(1:100,49)/AVOG,prof.plays(1:100,49),'+-');       grid on; xlabel('WV amt'); ylabel('plays [mb]'); legend('US STD','new prof','location','best')
  grid on;set(gca,'YDir','Reverse'); ylim([5e-3 1100])

figure(8);clf;semilogy(A.data(:,8),A.data(:,4)*STDATM,'o-',prof.gas_3(1:100,49)/AVOG,prof.plays(1:100,49),'+-');       grid on; xlabel('O3 amt'); ylabel('plays [mb]'); legend('US STD','new prof','location','best')
  grid on;set(gca,'YDir','Reverse'); ylim([5e-3 1100])

figure(9);clf;semilogy(A.data(:,5),A.data(:,4)*STDATM,'o-',prof.ptemp(1:100,49),prof.plays(1:100,49),'+-');       grid on; xlabel('T(z)'); ylabel('plays [mb]'); legend('US STD','new prof','location','best')
  grid on;set(gca,'YDir','Reverse'); ylim([5e-3 1100])

% -----------------------------------------------------------------
% write a new reference profile for use w/ SARTA in text form
if exist([sav.dir sav.fn_txt])
  fprintf(1,'oh oh .txt file already exists \n',[sav.dir sav.fn_txt])
  error('oh oh .txt file already exists')
end

FH = fopen([sav.dir sav.fn_txt],'w');
for j=1:length(A.textdata)
  fprintf(FH,'%s\n', A.textdata{j}');
end
%%% ! no altitude(m)  thickness(m) pressure(atm) temp(K) fixed=CO2     H2O        O3         CO         CH4        SO2       HNO3         N2O       NH3
for i=1:100; 
  ilay = 101-i;
  fprintf(FH,['%4d %11.5E %11.5E %11.5E %7.3f %11.5E %11.5E %11.5E %11.5E %11.5E' ...
     ' %11.5E %11.5E %11.5E %11.5E\n'], ...
     i, prof.palts(ilay,49), thick(ilay), prof.plays(ilay,49)/1013.25, prof.ptemp(ilay,49), ...
     prof.gas_2(ilay,49)/AVOG, prof.gas_1(ilay,49)/AVOG, prof.gas_3(ilay,49)/AVOG, ...
     prof.gas_5(ilay,49)/AVOG, prof.gas_6(ilay,49)/AVOG, prof.gas_9(ilay,49)/AVOG, ...
     prof.gas_12(ilay,49)/AVOG,prof.gas_4(ilay,49)/AVOG, prof.gas_4(ilay,49)/AVOG);
end
fclose(FH);
% --------------------------------------------------------------------
% For water continuum regression: write same output but in matlab format - 
% and ensure altitude order reversed:

for i=1:100; 
  ilay = 101-i;
  lev(ilay) = ilay;
  zref(ilay) = prof.palts(ilay,49);
  dzref(ilay) = thick(ilay);
  pref(ilay) = prof.plays(ilay,49);
  tref(ilay) = prof.ptemp(ilay,49);
  fref(ilay) = prof.gas_2(ilay,49)/AVOG;
  wref(ilay) = prof.gas_1(ilay,49)/AVOG;
  oref(ilay) = prof.gas_3(ilay,49)/AVOG;
  cref(ilay) = prof.gas_5(ilay,49)/AVOG;
  mref(ilay) = prof.gas_6(ilay,49)/AVOG;
  sref(ilay) = prof.gas_9(ilay,49)/AVOG;
  href(ilay) = prof.gas_12(ilay,49)/AVOG;
  nref(ilay) = prof.gas_4(ilay,49)/AVOG;
  aref(ilay) = prof.gas_11(ilay,49)/AVOG;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
FH = fopen('/home/chepplew/projects/sarta/cris_hr/refprof_regr49_1100_400ppm','r');
for i=1:9; info{i}=fgetl(FH); end
  newprof=textscan(FH,'%d %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(FH);
lev   = newprof{1}
zref  = newprof{2}(end:-1:1);
dzref = newprof{3}(end:-1:1);
pref  = newprof{4}(end:-1:1);
tref  = newprof{5}(end:-1:1);
fref  = newprof{6}(end:-1:1);
wref  = newprof{7}(end:-1:1);
oref  = newprof{8}(end:-1:1);
cref  = newprof{9}(end:-1:1);
mref  = newprof{10}(end:-1:1);
sref  = newprof{11}(end:-1:1);
href  = newprof{12}(end:-1:1);
nref  = newprof{13}(end:-1:1);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ! no altitude(m)  thickness(m) pressure(atm) temp(K) fixed=CO2     H2O        O3         CO         CH4        SO2       HNO3         N2O       NH3
savVars = {'zref','dzref','pref','tref','fref','wref','oref','cref','mref',...
           'sref','href','nref','aref'};
if exist([sav.dir sav.fn_mat])
  fprintf(1,'oh oh .mat file already exists \n',[sav.dir sav.fn_mat])
  error('oh oh .mat file already exists')
end
	   
save([sav.dir sav.fn_mat], savVars{:});

%%%%%%%%%%%%%%%% END of Reference Profile Work %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
sav.dir = '/home/sergio/git/ftc_dev/Refs/';
anew = load([sav.dir 'refprof_400ppm_aircraft_12km.mat']);
a0   = load([sav.dir 'refprof_regr49_1100_400ppm_b.mat']);
anew
a0
semilogy(a0.tref,a0.pref,'bo-',anew.tref,anew.pref,'r.-'); set(gca,'ydir','reverse'); ylim([min(a0.pref) max(a0.pref)]); ylabel('P [mb]');   xlabel('T(z)')
semilogy(a0.dzref,a0.pref,'bo-',anew.dzref,anew.pref,'r.-'); set(gca,'ydir','reverse'); ylim([min(a0.pref) max(a0.pref)]); ylabel('P [mb]'); xlabel('dz(z)')
loglog(a0.wref,a0.pref,'bo-',anew.wref,anew.pref,'r.-'); set(gca,'ydir','reverse'); ylim([min(a0.pref) max(a0.pref)]); ylabel('P [mb]');     xlabel('WV(z)')
loglog(a0.fref,a0.pref,'bo-',anew.fref,anew.pref,'r.-'); set(gca,'ydir','reverse'); ylim([min(a0.pref) max(a0.pref)]); ylabel('P [mb]');     xlabel('CO2(z)')
loglog(a0.oref,a0.pref,'bo-',anew.oref,anew.pref,'r.-'); set(gca,'ydir','reverse'); ylim([min(a0.pref) max(a0.pref)]); ylabel('P [mb]');     xlabel('O3(z)')
loglog(a0.mref,a0.pref,'bo-',anew.mref,anew.pref,'r.-'); set(gca,'ydir','reverse'); ylim([min(a0.pref) max(a0.pref)]); ylabel('P [mb]');     xlabel('CH4(z)')

%}
