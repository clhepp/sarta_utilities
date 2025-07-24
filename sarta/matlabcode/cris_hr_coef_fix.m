% function cris_hr_coef_fix.m
%
% purpose: manipulate Scott's CrIS hires fast code Coefficient sets to match
%     the actual CrIS spectral channels with 2 guard channels per edge.
%
% Source: /asl/data/sarta_database/Data_cris_nov04/Coef/
%
%
cd /home/chepplew/projects/sarta/cris_hr/

addpath /home/chepplew/projects/sarta/airs/code
addpath /asl/packages/ccast/source                   % seq_match

% ----------------------
% CrIS binary data files
% ----------------------
crspath = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/';
fnams   = {'set1.dat','set2.dat','set3.dat','set4_08.dat','set5_08.dat',...
           'set6_08.dat','set7_08.dat'};
%dfils = dir(strcat(crspath,'set*.dat'));
% 
clear crschns crsfrq crscoef;
lmerged = 0;

for mset = 1:7
  %fname   = strcat(crspath, sprintf('set%d.dat', mset));
  fname=strcat(crspath,fnams{mset});
  % read in the coefficients
  disp([num2str(mset) '  ' num2str(lmerged) '  ' fname]);
  [ichan, fchan, coef, info] = rdcoef( mset, lmerged,  fname);
    %whos ichan fchan coef info; disp(info);
  crsfrq{mset}  = fchan;
  crschns{mset} = ichan;
  crscoef{mset} = coef;
end
  
% record dimensions of the coef array.
[nfrq nlay ncoe] = size(crscoef{1});

% check total number of channels
sumcrs = 0;
for i = 1:7
  sumcrs = sumcrs + size(crscoef{i},1);
end

figure(1);clf;plot(crsfrq{1}, crscoef{1}(:,90,1),'.');grid on;title('CrIS');
figure(1);clf;hold on; grid on;
  for i=1:7 plot(crschns{i}, (i-1) + ones(1,numel(crschns{i})),'.'); end           
  legend('1','2','3','4','5','6','7','Location','northWest');
  



% Trim Scott's CrIS hi-res coefficient set to fit into 2223 (2 guard channel) freq grid.
% establish nominal HR grid:
band1 = [648.750:0.625:1096.250];
band2 = [1208.750:0.625:1751.250];
band3 = [2153.750:0.625:2551.250];
fchr  = [band1'; band2'; band3'];

% trim set 1: to match start freq. and sets 5 and 7 to match end freq.
trim1 = crsfrq{1}(23:end);                            % nchns: 491
trim5 = crsfrq{5}(1:73);                              % nchns:  73
trim7 = crsfrq{7}(1:25);                              % nchns:  25
% total number of channels now: 2215 (8 short) so find which channels are missing.
crhrfr = [trim1; crsfrq{2}; crsfrq{3}; crsfrq{4}; trim5; crsfrq{6}; trim7];

[Cx iC] = setdiff(fchr, crhrfr);          % <- 8 channels found
%   ch8= [1095.625 1096.250 1208.750  1209.375  1750.625  1751.250  2153.750  2154.375];
% add to: set2     set2     set3      set3      set3      set3      set4      set4

% -----------------------------------------------
% fabricate coefficients for the missing channels
clear extfrq ecoef;
msets = [2 3 4];
for m = msets % disp(m); end
  [nfrq nlay ncoe] = size(crscoef{m});
  if (m == 2) extfrq{m} = [crsfrq{2}; 1095.625; 1096.250]; end
  if (m == 3) extfrq{m} = sort([crsfrq{3}; 1208.750; 1209.375; 1750.625; 1751.250]); end
  if (m == 4) extfrq{m} = sort([crsfrq{4}; 2153.750; 2154.375]); end
  for k = 1:nlay
    for j = 1:ncoe
      ecoef{m}(:,k,j) = interp1(crsfrq{m},crscoef{m}(:,k,j),extfrq{m},'linear','extrap');
    end
  end
  fprintf('.');
end

% establish new ichans to match modified sets (retain original absolute channel mapping)
% and complete the frequencies of remaining sets.
clear echns;
[MX MY]  = seq_match(trim1, fchr);
echns{1} = MY;    extfrq{1} = fchr(MY);  ecoef{1} = crscoef{1}(23:end,:,:);  clear MX MY;
[MX MY]  = seq_match(extfrq{2}, fchr);
echns{2} = MY;    clear MX MY;
[MX MY]  = seq_match(extfrq{3}, fchr);
echns{3} = MY;    clear MX MY;
[MX MY]  = seq_match(extfrq{4}, fchr);
echns{4} = MY;    clear MX MY;
[MX MY]  = seq_match(trim5, fchr);
echns{5} = MY;    extfrq{5} = fchr(MY);  ecoef{5} = crscoef{5}(1:73,:,:);  clear MX MY;
[MX MY]  = seq_match(crsfrq{6}, fchr);
echns{6} = MY;    extfrq{6} = fchr(MY);  ecoef{6} = crscoef{6};            clear MX MY;
[MX MY]  = seq_match(trim7, fchr);
echns{7} = MY;    extfrq{7} = fchr(MY);  ecoef{7} = crscoef{7}(1:25,:,:);  clear MX MY;

% -----------------------------------------------
% write the new modified coefficient sets out
% -----------------------------------------------
dout = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/';
for m = 1:7
  fname = strcat(crspath,fnams{m});
  % read back in to get the info structure
  [ichan, fchan, coef, info] = rdcoef( m, lmerged,  fname);
  % modify parameters
  info.nchan = numel(echns{m});
  fout = sprintf('set%d_hrg2.dat',m);
  disp([num2str(m) '  ' num2str(lmerged) '  ' fnams{m} ' ' fout ' ' num2str(info.nchan)]);
  % 
  [iok] = wrtcoef( echns{m}, extfrq{m}, ecoef{m}, info, strcat(dout,fout) );

end

% variable CO2
fname = strcat(crspath,'co2_08.dat');
[ichan, fchan, coef, info] = rdcoef( 8, lmerged,  fname);       % <- enter 2 for CO2.
[MX MY]    = seq_match(fchan, fchr);
extfrq{8}  = fchan(MX);
ecoef{8}   = coef(MX,:,:);
info.nchan = numel(extfrq{8});
echns{8}   = MY;
%
fout = 'co2_hrg2.dat';
disp([num2str(8) '  ' num2str(lmerged) '  ' fname ' ' fout ' ' num2str(info.nchan)]);
% 
[iok] = wrtcoef( echns{8}, extfrq{8}, ecoef{8}, info, strcat(dout,fout) );

% optran
fname = strcat(crspath,'optran.dat');
[ichan, fchan, coef, info] = rdcoef( 9, 1,  fname);       % lmerged = 1
  whos ichan fchan coef info
[MX MY]    = seq_match(fchan, fchr);                      % appears to have 100% overlap.
% re-align channel IDs
echns{9} = MY;
fout = 'optran_hrg2.dat';
[iok] = wrtcoef( echns{9}, fchan, coef, info, strcat(dout,fout) );


% therm08.dat
fname = strcat(crspath, 'therm08.dat');
[ichan, fchan, coef, info] = rdcoef(11, 0,  fname);       % <- enter 0 for therm.
[MX MY]    = seq_match(fchan, fchr);                      % 2215 chans
extfrq{11} = sort([fchan(MX); ch8']);                     % add the missing 8.
echns{11}  = [1:2223]';
for k = 1:info.nlay
  for j = 1:info.ncoef
    ecoef{11}(:,k,j) = interp1(fchan,coef(:,k,j),extfrq{11},'linear','extrap');
  end
end
info.nchan = numel(extfrq{11});
fout = 'therm_hrg2.dat';
[iok] = wrtcoef( echns{11}, extfrq{11}, ecoef{11}, info, strcat(dout,fout) );

% ----------------------------------------------------------------------------
% create a tunmlt file (9 cols x 2223 rows. col1: echns, col2: fchr. then ones).

fout = 'tunmlt_ones2223.txt';
oFH  = fopen(strcat(dout,fout),'w');
for k = 1:2223
  fprintf(oFH,'  %d\t%8.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',...
     k,fchr(k),1,1,1,1,1,1,1);
end
fclose(oFH);

% write a subset copy of solar_data_cris08_08_08cm.txt
fin = 'solar_data_cris08_08_08cm.txt';
iFH = fopen(['/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Solar/' fin],'r');
Arr = fscanf(iFH,'%e %e %e', [3, Inf]);
fclose(iFH);
[MX MY] = seq_match(Arr(2,:),fchr);              % <- return 2215 values (8 short as above)
junk    = Arr(3,MX);
solar   = interp1(Arr(2,MX), junk, fchr);        % now have 2223 values
fout = 'solar_data_cris_hrg2.txt';
oFH = fopen(['/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Solar/' fout],'w');
for k = 1:2223
  fprintf(oFH, '  %d\t%8.3f\t%10.4e\n', k,fchr(k),solar(k));
end
fclose(oFH);

% figure(10);clf;plot(Arr(2,:),Arr(3,:),'.',fchr,solar,'.');grid on;

%}
