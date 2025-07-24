% compare_kc_od_spectroscopy
%
% compare the kCARTA Optical Depth breakout spectroscopy between different
%   HITRAN versions
% Takes original KCARTA predicts and does simple side-by-side comparison
%
% Note: arrays start at the SFC =layer 1

% Choose BreakOut
brkout = 'ozone';

% Choose original kCARTA OD paths
clear kc
kc.home = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
kc(1).src = [kc.home 'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];
kc(1).src = [kc.home 'REGR49_400ppm_H2016_May2021_AIRS2834_3CrIS_IASI/'];
kc(2).src = [kc.home 'REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/'];

% Choose target gas/breakout and assign breakout paths (bop)
% Ozone: from set2 (FOWP) b/o: F and FO, FWO, FWOP
bop.f      = 'F/';
bop.fo     = 'FO/';
bop.fwo    = 'FWO/';
bop.fwop   = 'FWOP/';

bop.ch4f   = 'ch4bandF/';       % == bop.f (F/)
bop.ch4    = 'ch4bandCH4/';
bop.so2    = 'so2bandS/';

% and get individual profile file listings (nominally 49 each)
clear flist;
flist(1).f    = dir([kc(1).src bop.f '*.mat']);
flist(1).fo   = dir([kc(1).src bop.fo '*.mat']);
flist(1).fwo  = dir([kc(1).src bop.fwo '*.mat']);
flist(1).fwop = dir([kc(1).src bop.fwop '*.mat']);
flist(1).ch4  = dir([kc(1).src bop.ch4 '*.mat']);
flist(1).so2  = dir([kc(1).src bop.so2 '*.mat']);

flist(2).f    = dir([kc(2).src bop.f '*.mat']);
flist(2).fo   = dir([kc(2).src bop.fo '*.mat']);
flist(2).fwo  = dir([kc(2).src bop.fwo '*.mat']);
flist(2).fwop = dir([kc(2).src bop.fwop '*.mat']);
flist(2).ch4  = dir([kc(2).src bop.ch4 '*.mat']);
flist(2).so2  = dir([kc(2).src bop.so2 '*.mat']);

% Load values into memory
clear vals
ipr = 1;
vals(1).f    = [];
vals(1).fo   = [];
vals(1).fwo  = [];
vals(1).ch4  = [];
vals(1).so2  = [];
for i=1:length(flist(1).f)
  junk = load([flist(1).f(i).folder '/' flist(1).f(i).name],'riasi_all');
  vals(1).f(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(1).fwo(i).folder '/' flist(1).fwo(i).name],'riasi_all');
  vals(1).fwo(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(1).fo(i).folder '/' flist(1).fo(i).name],'riasi_all');
  vals(1).fo(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(1).ch4(i).folder '/' flist(1).ch4(i).name],'riasi_all');
  vals(1).ch4(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(1).so2(i).folder '/' flist(1).so2(i).name],'riasi_all');
  vals(1).so2(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
%
  ipr=ipr+14;
  fprintf(1,'.')
end

ipr = 1;
vals(2).f    = [];
vals(2).fo   = [];
vals(2).fwo  = [];
vals(2).ch4  = [];
vals(2).so2  = [];
for i=1:length(flist(2).f)
  junk = load([flist(2).f(i).folder '/' flist(2).f(i).name],'riasi_all');
  vals(2).f(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(2).fwo(i).folder '/' flist(2).fwo(i).name],'riasi_all');
  vals(2).fwo(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(2).fo(i).folder '/' flist(2).fo(i).name],'riasi_all');
  vals(2).fo(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(2).ch4(i).folder '/' flist(2).ch4(i).name],'riasi_all');
  vals(2).ch4(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
  junk = load([flist(2).so2(i).folder '/' flist(2).so2(i).name],'riasi_all');
  vals(2).so2(:,:,ipr:ipr+13) = permute(junk.riasi_all,[3 2 1]);
%
  ipr=ipr+14;
  fprintf(1,'.')
end

% get the frequency grid
load([flist(2).f(i).folder '/' flist(2).f(i).name],'fiasi')

%{
 phome = '/home/chepplew/figs_sync/';
 diff1 = squeeze(vals(1).f(90,:,1)) - squeeze(vals(1).fo(90,:,1));
 diff2 = squeeze(vals(2).f(90,:,1)) - squeeze(vals(2).fo(90,:,1));
 diff1 = squeeze(vals(1).f(9,:,1)) - squeeze(vals(1).ch4(9,:,1));
 diff2 = squeeze(vals(2).f(9,:,1)) - squeeze(vals(2).ch4(9,:,1));
 dfrc12 = (diff1 - diff2)./diff1;
 fh1=figure('visible','off');
 plot(fiasi, diff1,'-', fiasi, diff2,'-');grid on; legend('diff1','diff2');
 plot(fiasi, dfrc12,'-');grid on;xlabel('wvn cm-1');ylabel('fractional trans diff')
    title('Ozone lev:90. prof:1');axis([650 2800 -10 10]);
 saveas(fh1, [phome 'plot_20220712_00c.png'],'png')
%}

