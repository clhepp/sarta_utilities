function [iok] = write_tunmlt_ones(csens)

% create a tunmlt file (9 cols x 2223 rows. col1: echns, col2: fchr. then ones).


% Check csens
csens = upper(csens);
if(~ismember(csens,['CRIS_LR','CRIS_HR','CHIRP','AIRS_L1C','IASI']))
  error('Invalid sensor type')
  return
end

% default destination
dout = '/home/chepplew/data/sarta/';

% Default header
clear hdr;
hdr{1} = ['% ' csens];
hdr{2} = '%';
hdr{3} = ['% no tuning'];
hdr{4} = '%';
hdr{5} = ['%  ID  freq          fixed   w lines w con   ozone    CO      CH4    nonLTE'];
hdr{6} = ['%---- --------       ------  ------- ------  ------- ------- ------- -------'];


switch csens
  case 'CRIS_LR'
    nchan = 1329;
    fnout = 'tunmlt_ones_cris_lr_g4.txt';
    load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    % check vchan and idchan
    % TBD
    if(length(vchan) ~= nchan | length(idchan) ~= nchan)
      error('Unexpected spectral grid length')
      return;
    end

  case 'CRIS_HR'
    nchan = 2235;
    fnout = 'tunmlt_ones_cris_hr_g4.txt';

  case 'CHIRP'
    nchan = 1609;
    fnout = 'tunmlt_ones_chirp.txt';

  case 'AIRS_L1C'
    nchan = 2645;
    fnout = 'tunmlt_ones_airs_l1c.txt';

  case 'IASI'
    nchan = 8461;
    fnout = 'tunmlt_ones_iasi.txt';

end

% Write to file
oFH  = fopen(strcat(dout,fnout),'w');
for k = 1:length(hdr)
  fprintf(oFH, '%s\n', hdr{k});
end
for k = 1:nchan
  fprintf(oFH,'  %d\t%8.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',...
     idchan(k),vchan(k),1,1,1,1,1,1,1);
end
iok = fclose(oFH);

% iok = 0;
