% Match up CrIS thermal FTC coefficients

addpath /home/chepplew/projects/sarta/matlabcode

% ------------------------
% load Scott's coefficient
% ------------------------
FNCOF11 = '/asl/data/sarta_database/Data_cris_nov04/Coef/therm08.dat';
[ich11 fch11 coef11 info11] = rdcoef(11,0,FNCOF11);

% load CrIS grid

xx=load('/home/chepplew/gitLib/ftc_dev/chanLists/list_cris_hrg4');
ichan = xx(:,1);
freq  = xx(:,2);
clear xx;

% 

[IX IY] = intersect(fch11, freq);
[JX JY] = intersect(freq, fch11); 

therm_coef = coef11(IY,1,:);
therm_freq = fch11(IY);
therm_chan = JY;

therm_info = struct;
therm_info = info11;
therm_info.nchan = numel(therm_chan);

% write fortran binary coefficient file

fname = 'therm_matched.dat';
[iok] = wrtcoef(therm_chan, therm_freq, therm_coef, therm_info, fname);

% ------------------------------------
% load Scott's CrIS hr CO2 coefficient:
% ------------------------------------
FNCOF8 = '/asl/data/sarta_database/Data_CrIS_apr09/Coef/co2g4.dat';
[ich8 fch8 coef8 info8] = rdcoef(8,0,FNCOF8);

[IX IY] = intersect(fch8, freq);
[JX JY] = intersect(freq, fch8); 
co2_coef = coef8(IY,:,:);
co2_freq = fch8(IY);
co2_chan = JY;

co2_info = struct;
co2_info = info8;
co2_info.nchan = numel(co2_chan);

% write fortran binary coef file:
fname = 'co2_matched.dat';
[iok] = wrtcoef(co2_chan, co2_freq, co2_coef, co2_info, fname);



