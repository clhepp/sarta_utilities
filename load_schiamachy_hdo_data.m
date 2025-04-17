

% HDO_IMPAv20 Monthly
d.home = '/home/chepplew/data/SCHIAMACHY/HDO_IMAPv20/2006/'; 
d.list = dir([d.home '/*/*/SCIA_IMAP_*_small.nc'])

fname = [d.list(1).folder '/' d.list(1).name];

fninfo = ncinfo(fname)

% Groups.Name={'OrbitGeometry','OrbitHeader','RetrievalResults','HighLevelResults'}

ngrps = length(fninfo.Groups);
for ngrp = 1:ngrps
  if(ismember(ngrp,[1 2])) 
    for i = 1:length(fninfo.Groups(ngrp).Variables)
      vars(ngrp).nams{i} = fninfo.Groups(ngrp).Variables(i).Name;
    end
  end
  if(ismember(ngrp,[3 4]))
    for i = 1:length(fninfo.Groups(ngrp).Groups.Variables)
      vars(ngrp).nams{i} = fninfo.Groups(ngrp).Groups.Variables(i).Name;
    end
  end
end

% [21 3 1 29]

lat = []; lon = []; tim = []; deltaD = []; hdo = []; ch4 = []; h2o = [];
xh2o = []; qflag = [];

for ifn=1:length(d.list)

  fname = [d.list(ifn).folder '/' d.list(ifn).name];

  lat   = [lat; ncread(fname,'/OrbitGeometry/center_latitude')];
  lon   = [lon; ncread(fname,'/OrbitGeometry/center_longitude')];
  tim   = [tim;  ncread(fname,'/OrbitGeometry/time_mjd')];

  deltaD = [deltaD;  ncread(fname,'/HighLevelResults/hdo_v20/deltaD')];
  h2o    = [h2o; ncread(fname,'/HighLevelResults/hdo_v20/h2o_vcd')];
  hdo    = [hdo; ncread(fname,'/HighLevelResults/hdo_v20/hdo_vcd')];
  ch4    = [ch4; ncread(fname,'/HighLevelResults/hdo_v20/ch4_vcd')];

  xh2o    = [xh2o; ncread(fname,'/HighLevelResults/hdo_v20/XH2O_retrieved')];
  
  qflag   = [qflag; ncread(fname,'/HighLevelResults/hdo_v20/QualityFlag')];
  
fprintf(1,'.')
end

whos lat lon tim deltaD ch4 h2o hdo xh2o

iiwnt = find(qflag ==1);
simplemap(lat(iiwnt), lon(iiwnt), xh2o(iiwnt) );
clf;simplemap(lat(iiwnt), lon(iiwnt), deltaD(iiwnt) )
  load llsmap5
  colormap(llsmap5); caxis)[-1000 1000])

% --------------------------------------------------------------
% HDO_SICORv10 orbit files (approx 414 per month)
% --------------------------------------------------------------
d.home = '/home/chepplew/data/SCHIAMACHY/HDO_SICORv10/2006/06/';
d.list = dir([d.home 'sci_l2_h2o_hdo*.nc']);

ifn=9;
fname = [d.list(ifn).folder '/' d.list(ifn).name]

fninfo = ncinfo(fname);
ngrps = length(fninfo.Groups)

%  %{'diagnostics','instrument','meteo','side_product','target_product'} 

clear vars
for ngrp = 1:ngrps
  for i = 1:length(fninfo.Groups(ngrp).Variables)
    vars(ngrp).nams{i} = fninfo.Groups(ngrp).Variables(i).Name;
  end
end

lat = []; lon = []; tim = []; spres = []; plevs = []; ptemp = [];
h2o_col = []; hdo_col = [];

for ifn = 9:23
  fname = [d.list(ifn).folder '/' d.list(ifn).name];

  lat   = [lat;   ncread(fname,'/instrument/latitude_center')];
  lon   = [lon;   ncread(fname,'/instrument/longitude_center')];
  tim   = [tim;   ncread(fname,'/instrument/time')'];
  spres = [spres; ncread(fname,'/meteo/surface_pressure')];
  plevs = [plevs; ncread(fname,'/meteo/pressure_levels')'];
  ptemp = [ptemp; ncread(fname,'/meteo/temperature')'];

  h2o_col = [h2o_col; ncread(fname,'/target_product/h2o_column')];
  hdo_col = [hdo_col; ncread(fname,'/target_product/hdo_column')];
 
  fprintf(1,'.')
end

whos lat lon tim spres plevs ptemp *_col

% deltaD = 1000*(VCD(HDO)/VCD(H2O) - 1) - 20  % mil

deltaD = 1000.0*((abs(hdo_col)./abs(h2o_col)) - 1.0);

