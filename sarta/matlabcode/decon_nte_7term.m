% decon_nte_7term

% script to translate the IASI non-LTE coefficient set from IASI to CrIS grid

% Note: the original IASI coefficient set spans only from 2215 to 2392 wn w/ 0.25 wn res.
% CrIS band 3 spans 2155 to 2550 wn (hi-res mode = 0.625 wn).

% This version to be run manually in matlab. NOT for general application.

% The write section is complete but commented.

% Draft v 1.  27-Apr-2016.
% Author: Hepplewhite

cd /home/chepplew/projects/sarta/cris_hr/

addpath /asl/packages/iasi_decon
addpath /asl/packages/ccast/source                        % inst_params.m
addpath /home/chepplew/projects/sarta/airs/code           % rdcoef.m

% defaults
opt1.hapod   = 0;         % no Hamming apodization (my preference - not much difference)
opt1.nguard  = 0;         % no guard channels
opt1.resmode = 'hires2';  % 0.625 wn all 3 bands.

% load hi-res CrIS grid w/ 2 guard channels
xx = load('fcris_hrg2.mat'); fcris = xx.fcris_hrg2; clear xx;

% load original IASI data
fn = '/asl/data/sarta_database/Data_IASI_sep08/Coef/nte_7term.dat';
[ichan, fchan, coef, info] = rdcoef(11, 0, fn);     % respond: 1, 7.

%{
figure(1);clf;plot(fchan,coef(:,1,2),'.-');
%}

% prep deconvolve to CrIS
% IASI params
iasi.v1 = 645;            % iasi band low
iasi.v2 = 2760;           % iasi band high
iasi.dv = 0.25;           % IASI dv

% CrIS params
wlaser = 773.1301;        % nominal value
bstr{1} = 'LW';           % band by number
bstr{2} = 'MW'; 
bstr{3} = 'SW';

% initialize outputs
cdec = []; fdec = []; 
bpts = []; bv1 = []; bv2 = [];

% CrIS band 3 (SW)
bi = 3;

band = bstr{bi};
[inst, user] = inst_params(band, wlaser, opt1);

% get the passband plus usable filter wings
tv1 = max(iasi.v1, user.v1 - user.vr);
tv2 = min(iasi.v2, user.v2 + user.vr);
tvr = min(user.v1 - tv1, tv2 - user.v2);

% trim IASI data to the filter span 
ix = find(tv1 <= fchan & fchan <= tv2);
ftmp = fchan(ix); 

[m, n, p]   = size(coef);           % expected: m=687, n=1, p=7.
for ip = 1:p
  ctmp = coef(ix, 1, ip); 

  % apply the bandpass filter
  ctmp = bandpass(ftmp, ctmp, user.v1, user.v2, tvr);

 % convolve to the CrIS user grid
  [r1, f1] = iasi_decon(ctmp, ftmp, user.dv, opt1);
  fdec = f1; ctmp = r1; clear r1;
  cdec(:,n,ip) = real(ctmp); clear ctmp;
end

% record common channels (for writing new coefficient file)
[si sj] = seq_match(fdec, fcris);

%{
%figure(2);clf;plot(fdec,cdec(:,1,2),'.-');grid on;
for ip = 1:7
  figure(ip);clf;plot(fchan,coef(:,1,ip),'.-',fdec,cdec(:,1,ip),'r.-');grid on;
    legend('iasi','i2c');title(['nte.7term i2c coef' sprintf('%d',ip)]);
  %saveas(gcf,['./figs/nte_7term_i2c_coef' sprintf('%d',ip) '.png'],'png');
  pause(5);
end
%}

% option to write new coefficient set out:
dp = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/';
fo = 'nte_7term_hrg2.dat';

o_info = struct;
o_info.set       = info.set;
o_info.ncoef     = info.ncoef;
o_info.nlay      = info.nlay;
o_info.gasid     = info.gasid;
o_info.nbreakout = info.nbreakout;
o_info.ncoefeach = info.ncoefeach;
o_info.startin   = info.startind;
o_info.nchan     = numel(fdec);

%{
fid = fopen(strcat(dp, fo),'w','ieee-be');

% Value of FORTRAN record marker for each channel
ifm = round( 4*(1 + 1 + o_info.nlay*o_info.ncoef) ); % exact integer

% Loop over the channels
for k = 1:o_info.nchan

   % FORTRAN start-of-record marker
   icount = fwrite(fid,ifm,'integer*4');

   % data for this channel
   icount = fwrite(fid,sj(k),'integer*4');     % index/channel number for 2223 hr_g2 set.
   icount = fwrite(fid,fdec(k),'real*4');
   for i = 1:o_info.nlay
      icount = fwrite(fid,cdec(k,i,:),'real*4');
   end

   % FORTRAN end-of-record marker
   icount=fwrite(fid,ifm,'integer*4');

end % end of loop over bands

fclose(fid);
%}

