function [iok] = write_solardata_for_sarta(csens)

% from 605 to 2430 cm-1 use /asl/rta/kcarta/solar/srad*.mat
% from 2430 to 3000 cm-1 use atmos 
%      /asl/s2/hannon/Solar_data/ATMOS/atmos_corrected_solar.mat
%    fortunately both spectra are on same 0.0025 cm-1 grid.
%
% CLH v2: use: /asl/s2/hannon/Solar_data/Chris_Barnet/fine_solar.mat
%      fout, rout 

% Check csens
csens = upper(csens);
if(~ismember(csens,['CRIS_LR','CRIS_HR','CHIRP','AIRS_L1C','IASI']))
  error('Invalid sensor type')
  return
end

% default destination
dout = '/home/chepplew/data/sarta/';

% default header
hdr = [];
hdr{1} = ['!chan frequency  solar'];
hdr{2} = ['!----  -------  -------'];

% Get original solar data
%{
dfreq  = 0.0025;    % cm-1
d.home = '/asl/rta/kcarta/solar/';
d.dir  = dir([d.home 'srad*.mat']);

srad = [];
for i = 1:length(d.dir)
  x    = load([d.dir(i).folder '/' d.dir(i).name]);
  srad = [srad x.srad];
end
sfreq = [605.000:0.0025:2829.9975];
% truncate to 2430 cm-1
iix = find(sfreq>2430,1)-1;

atmos = load('/asl/s2/hannon/Solar_data/ATMOS/atmos_corrected_solar.mat');
iiy = find(atmos.fsol>2430,1);

xrad = [srad(1:iix) atmos.rsol(iiy:end)];

%}

% fout, rout
load('/asl/s2/hannon/Solar_data/Chris_Barnet/fine_solar.mat');



switch csens
  case 'CRIS_LR'
    nchan = 1329;
    fnout = 'solardata_cris_lr_g4.txt';
    load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    % check vchan and idchan
    % TBD
    if(length(vchan) ~= nchan | length(idchan) ~= nchan)
      error('Unexpected spectral grid length')
      return;
    end  
  
  case 'CRIS_HR'
  
  case 'CHIRP'
  
  case 'AIRS_L1C'
  
  case 'IASI'
  
end

% Interpolate original data to the sensor grid
%%% irad = interp1(sfreq, srad, vchan,'linear');
irad = interp1(fout, rout, vchan,'linear');

% Write to file
oFH  = fopen(strcat(dout,fnout),'w');
for k = 1:length(hdr)
  fprintf(oFH, '%s\n', hdr{k});
end
for k = 1:nchan
  fprintf(oFH,'  %d\t%8.3f\t%6.3f\n',...
     idchan(k),vchan(k),irad(k));
end
iok = fclose(oFH);

% iok = 0;
 
%{
plot(sfreq, srad, '-')

%}
