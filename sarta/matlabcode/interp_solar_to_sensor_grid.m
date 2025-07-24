% function [] = interp_solar_to_sensor_grid(csens)

% INPUT: csens strig {'cris_fsr','iasi','airs_l1c','chirp'}
% OUTPUT

% Source solar spectrum data
fn_src = '/home/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Solar/solar_data_cris08_08_08cm.txt';
fn_src = '/asl/s2/hannon/Solar_data/Chris_Barnet/nast_solar_data.txt';
  %   column: 1:line, 2:wvn, 3:radiance. Spectral range: 605.0 to 2830.0 cm-1.  rad units W.m-2.sr-1.cm

% Original Scott header 
fn_hdr = '/asl/s2/hannon/AIRS_prod08/Fit_ftc/Solar/header';

% Check input argument
csens = upper(csens)
if(~ismember(csens,{'AIRS_L1C','CRIS_FSR','IASI','CHIRP'})) error('Invalid sensor'); return; end

% Get sensor grid
switch csens
  case 'CHIRP'
    fn_grid = '/home/chepplew/gitLib/ftc_dev/chanLists/chirp_1702_list_all';
    fn_out  = '/home/chepplew/data/sarta/prod_2020/chirp/feb2020/dbase/Solar/solardata.txt';
  case 'IASI'
    fn_grid = '';
  case 'AIRS_L1C'
    fn_grid = '';
  case 'CRIS_FSR'
    fn_grid = '';
end

% Load grid data
dgrid = importdata(fn_grid);

% load source solar spectrum
dsrc = importdata(fn_src,' ',5);
[src_wvn, iz] = unique(dsrc.data(:,2));
src_rad       = dsrc.data(iz,3);

 
% Interpolate to desired grid
intsol   = interp1(src_wvn, src_rad, dgrid(:,2));        % now have 2223 values

% Write data to file
% check output directory exists
dr_out = fileparts(fn_out);
if(~exist(dr_out)) error('Output directory does nto exist - plz check'); return; end

FH = fopen(fn_out,'w')

  for i = 1:length(dsrc.textdata)
    fprintf(FH, '%s\n', str2mat(dsrc.textdata{i})); 
  end
%
  for i=1:length(intsol)
    fprintf(FH,'  %6d  %8.4f  %8.3f\n', i, dgrid(i,2), intsol(i));
  end

fclose(FH);

%{
figure(1);clf;plot(dsrc(:,2), dsrc(:,3),'.-')
  hold on; plot(dgrid(:,2), intsol,'.-')
%}
