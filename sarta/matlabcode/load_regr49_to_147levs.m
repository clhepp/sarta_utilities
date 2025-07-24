% function [prof49] = load_regr49_to_147levs(opts)
%
% Load the raw regr49 levels data, interpolate to 
%      the 147 levels of the SAF profiles
%      for later concatenation
% update with AFGL and SABER atomicO
%
% INPUT:
%       opts structure with fields:
%





% ---------------------------------------------------------------
% for splicing in missing upper profiles use afgl, load in now
% ---------------------------------------------------------------
afgl.a1 = quick_read_afgl(1,1);
afgl.a2 = quick_read_afgl(1,2);
afgl.a3 = quick_read_afgl(1,3);
afgl.a4 = quick_read_afgl(1,4);
afgl.a5 = quick_read_afgl(1,5);
afgl.a6 = quick_read_afgl(1,6);
afgl.tstd = [afgl.a1.tstd, afgl.a2.tstd, afgl.a3.tstd, afgl.a4.tstd, ...
             afgl.a5.tstd, afgl.a6.tstd];

% --------------------------------------------------------
% add regr49 reference profiles - interpolated to plevs147
% NB:  plevs147 TOA->SFC, prof49.plevs SFC->TOA => reverse 
% --------------------------------------------------------
fn.r49_op = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
           'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];
fn.r49_ip = ['/asl/s1/sergio/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
           'REGR49_H2012_June2016/regr49_1100.ip.rtp'];
% !! check head.gunit :  g/g 21. !!
% NB Only the first 50 levels has values !!
[head49,~,prof49,~] = rtpread(fn.r49_ip);

% load RMM for gases
rmms = importdata('/home/chepplew/myLib/data/RMM_of_HITRAN_gases.txt');

% convert gas units from 10 (mmr) to 21) ppmv
for ii = 1:length(head49.glist)
  gasID = head49.glist(ii);
  gasUI = head49.gunit(ii);
  gasRM = rmms.data(gasID);
  gasAM = prof49.(['gas_' num2str(gasID)]);
  y     = toppmv(prof49.plevs, prof49.ptemp, gasAM, gasRM, gasUI);
  prof49.(['gas_' num2str(gasID)]) = y;
end

% update co2ppm to 400 ppm (mmr 611.11 g/g dry air)
% co2_fac = 6.111e-4/max(prof49.gas_2(:));      % 6.111/5.62   (*10-4  g/g)
co2_fac = 400.0/max(prof49.gas_2(:));      % ppmv
max(prof49.gas_2(:))
prof49.co2ppm = 400.0*ones(size(prof49.co2ppm));
prof49.gas_2  = co2_fac * prof49.gas_2;

% Set plon to 160.0 deg, get salti and landfrac then force to < 100 m
prof49.plon = 160.0*ones(1,49);
[salti, landfrac] = usgs_deg10_dem(prof49.plat, prof49.plon);
iibig = find(salti > 100.0)
if(~isempty(iibig))
  salti(iibig) = 100.0;
end
prof49.salti    = salti;
prof49.landfrac = landfrac;
% set spres = 1013.25 hPa.
prof49.spres    = 1013.25*ones(1,49);

% -------------------------------------------------
% prof49 profiles are mixture of lengths and ranges
%        get min plevs and splice on afgl if needed
% -------------------------------------------------
clear p49min i49min
prof49.plevs(prof49.plevs == 0) = NaN;

% for short profiles extend using afgl (NB profiles are SFC->TOA)
for ip = 1:49
  [p49min(ip), i49min(ip)] = min(prof49.plevs(:,ip));
  if(p49min(ip) > 4.0E-5)
% find best matched AFGL zone to use assigned to IZ
    iipA = find(afgl.a1.pstd <= p49min(ip));
    %%iipB = find(prof49.plevs(:,ip) < p49min(ip),1);
    %%iipB = find(prof49.plevs(:,ip) > p49min(ip),1,'last');
    iipB = i49min(ip)-1;         % ! -1

    [~, IZ] = min( abs(afgl.tstd(iipA(1),:) - prof49.ptemp(i49min(ip),ip) ));
    fprintf(1,'ip, pmin(ip): iipB,IZ %4i, %5.2e, %4i,%4i\n', ...
            ip, p49min(ip), iipB,IZ)

    eval( ['Tjunk=afgl.a' num2str(IZ) '.tstd(:);'] );
    eval( ['Pjunk=afgl.a' num2str(IZ) '.pstd(:);'] );
    T_off = prof49.ptemp(iipB+1,ip) - Tjunk(iipA(1));
    prof49.ptemp(iipB+[1:length(iipA)],ip) = Tjunk(iipA) + T_off;
    prof49.ptemp(prof49.ptemp == 0) = NaN;

    afgl.g1 = quick_read_afgl(1,IZ);
    g1_off  = prof49.gas_1(iipB+1,ip) - afgl.g1.qiAtm(iipA(1));
    prof49.gas_1(iipB+[1:length(iipA)],ip) = afgl.g1.qiAtm(iipA) + g1_off;

    afgl.g2 = quick_read_afgl(2,IZ);
    g2_off  = prof49.gas_2(iipB+1,ip) - afgl.g2.qiAtm(iipA(1));
    prof49.gas_2(iipB+[1:length(iipA)],ip) = afgl.g2.qiAtm(iipA);

    afgl.g3 = quick_read_afgl(3,IZ);
    g3_off  = prof49.gas_3(iipB+1,ip) - afgl.g3.qiAtm(iipA(1));
    prof49.gas_3(iipB+[1:length(iipA)],ip) = afgl.g3.qiAtm(iipA);

    afgl.g5 = quick_read_afgl(5,IZ);
    prof49.gas_5(iipB+[1:length(iipA)],ip) = afgl.g5.qiAtm(iipA);

    afgl.g6 = quick_read_afgl(6,IZ);
    prof49.gas_6(iipB+[1:length(iipA)],ip) = afgl.g6.qiAtm(iipA);

    afgl.g9 = quick_read_afgl(9,IZ);
    prof49.gas_9(iipB+[1:length(iipA)],ip) = afgl.g9.qstd(iipA);
% if the field is extended all other field extended values are set to zero
%    so the test at the start of this loop will fail.
    prof49.plevs(iipB+[1:length(iipA)],ip) = Pjunk(iipA);
    prof49.plevs(prof49.plevs == 0) = NaN;   
    prof49.nlevs(ip)  = length(find(~isnan(prof49.plevs(:,ip))));
  end
  fprintf(1,'.')
end

% --------------------------------------------------
% interpolate to the 147 plevs and revserse TOA->SFC
% --------------------------------------------------
%%plevs147  from subset_saf() load(d.fin)
prof       = prof49;
prof.ptemp = []; prof.gas_1 = []; prof.gas_2 = [];
prof.gas_3 = []; prof.gas_5 = []; prof.gas_6 = [];
prof.gas_9 = []; prof.plevs = []; prof.nlevs = [];
prof.plevs = repmat(plevs147,1,49);
for ip = 1:49
  nlev   = prof49.nlevs(ip);
  plevg  = log(prof49.plevs(1:nlev,ip));
  psurf  = prof49.plevs(1,ip);
  prof.ptemp(:,ip) = interp1(plevg, prof49.ptemp(1:nlev,ip), log(plevs147), ...
        [],'extrap');
  prof.gas_1(:,ip) = interp1(plevg, prof49.gas_1(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_2(:,ip) = interp1(plevg, prof49.gas_2(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_3(:,ip) = interp1(plevg, prof49.gas_3(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_5(:,ip) = interp1(plevg, prof49.gas_5(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_6(:,ip) = interp1(plevg, prof49.gas_6(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_9(:,ip) = interp1(plevg, prof49.gas_9(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
% check range limits. SFC pressure
  prof.gas_1(prof.gas_1 < 0) = 0;
  prof.gas_2(prof.gas_2 < 0) = 0;
  prof.gas_3(prof.gas_3 < 0) = 0;
  iiovr = find(plevs147 > prof49.plevs(1,ip));
  iirng = find(plevs147 <= prof49.plevs(1,ip));
  prof.nlevs(ip)   = length(iirng);
  %prof.plevs(:,ip) = plevs147(1:length(iirng));
  prof.plevs(iiovr,ip) = 0.0;
end

% ----------- assign needed fields to prof (49) ------
prof.rlat   = prof.plat;
prof.rlon   = prof.plon;
prof.rtime  = repmat(dtime2tai('2020/06/01'), 1,49);
prof.whichbunch = 7.0 * ones(1,49);

% ------------------------------------------
% add atomicO to prof49, are all in June
% ------------------------------------------
%  Load SABER data if not already loaded
% ---------------------------------------------------------------
% get afgl atomicO (gas.34)
afgl34 = quick_read_afgl(34,1);

if(~exist(atmonicO))
  % get Saber atomicO
  fn.saber = ['/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES_FOR_MANUEL/' ...
            'AtomicO_SABER_Mlynczak/atox_athy_night_YY2014_V1.02.nc'];
  [atomicO,~]  = read_netcdf_h5(fn.saber);
  atomicO.ppmv = toppmv(atomicO.pressure*ones(1,length(atomicO.lat)), ...
               atomicO.ktemp,atomicO.qatox,16,12);
  atomicO.lon = wrapTo180(atomicO.lon);

  % Match time/lat/lon
  daysINyear  = [0 31 28 31 30 31 30 31 31 30 31 30 31];
  cdaysINyear = cumsum(daysINyear);
  for ii = 1 : 12
    boo = find(atomicO.day > cdaysINyear(ii) & atomicO.day <= cdaysINyear(ii+1));
    atomicO.month(boo) = ii;
  end
end
% ----------------------------------------------------------
clear saber_*
mons = month(tai2dtime(prof.rtime));
ii = 6;
  booO = find(atomicO.month == ii);
  booP = find(mons == ii);
  [A,B,C,D,E] = haversine3_matrix(prof.rlat(booP), prof.rlon(booP), ...
                atomicO.lat(booO), atomicO.lon(booO), 5E7);
  for jj = 1:length(booP)
     [saber_O(:,booP(jj)), saber_p] = add_afgl_g34(atomicO.ppmv(:,booO(B(jj))), ...
           atomicO.pressure, afgl34);
  end
% interpolate to prof.plevs
for ii = 1 : length(prof.rlat)
  pjunk = prof.plevs(:,ii);
  saber_O_prof(:,ii) = interp1(log(saber_p),saber_O(:,ii),log(pjunk));
end
% replace existing zero values in prof.gas_34
prof.gas_34 = saber_O_prof;

