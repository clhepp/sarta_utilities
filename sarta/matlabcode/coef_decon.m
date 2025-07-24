function [coef info ccoef o_info] = coef_decon(cset, lmerged);


% function [ ] = coef_decon()
%
% INPUT:
%       cset:    integer [1 to 7] (was 'set' but this conflicts w/ Matlab libs!!!)
%       fname:  original coefficient binary data file (path hard wired)
%
% OUTPUT:
%       coef_orig:   original IASI coefficients
%       coef_dec:    deconvolved onto CrIS grid
%       info:        structure returnd from rdcoef.m
%
% PURPOSE:
%       translate Scott's IASI SARTA coefficients to the CrIS spectral grid
%
% ASSUMPTIONS:
%       For use with IASI coefficients sets 1 to 7 from 
%       /asl/data/sarta_database/Data_IASI_may09/Coef/
%
% HISTORY:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% work in my directory
cd /home/chepplew/projects/sarta/airs/code

% function and data paths 
addpath /asl/packages/iasi_decon/
addpath /asl/packages/ccast/source                             % inst_params
addpath /asl/packages/airs_decon/source                        % hamm_app

% check input params
if(cset > 12 || cset < 1) fprintf(1,'Error: cset not valid\n'); end

% default and other variables
%lmerged = 0;
load('/asl/data/iremis/danz/iasi_f.mat');                      % fiasi [8461x1]
xx = load('/home/chepplew/gitLib/asl_sno/data/cris_freq_2grd.mat');  
  fcris = xx.vchan; clear xx;                                  % 1317 chns (12 guard chans)

% expected binary data files
dpath = '/asl/data/sarta_database/Data_IASI_may09/Coef/';
dfils = dir(strcat(dpath,'set*.dat'));
othfn = {{8.1,'co2.dat'}, {8.2,'so2.dat'}, {8.3, 'hno3.dat'}, ...
         {9, 'optran.dat'}, {10, 'n2o.dat'}};
if(cset >=1 && cset <= 7)
   mset  = fix(cset);
   fname = strcat(dpath,dfils(mset).name);
end
if(cset >= 8 && cset <= 8.3)
  mset = 8;
  switch cset
    case 8.1
      k = 1; fname = strcat(dpath, cell2mat(othfn{k}(2)));
    case 8.2
      k = 2; fname = strcat(dpath, cell2mat(othfn{k}(2)));
    case 8.3
      k = 3; fname = strcat(dpath, cell2mat(othfn{k}(2)));
  end
end
if(cset == 9 || cset == 10)
  mset = fix(cset); k=4 + mset-9; fname = strcat(dpath, cell2mat(othfn{k}(2)));
end
  
% read in the coefficients
disp([num2str(mset) '  ' num2str(lmerged) '  ' fname]);
[ichan, fchan, coef, info] = rdcoef( mset, lmerged,  fname);
  whos ichan fchan coef info;
  disp(info);
  
% record dimensions of the coef array.
[nfrq nlay ncoe] = size(coef);

% prep for decon
opt1       = struct;
opt1.hapod = 0;

clear ccoef;
for ic = 1:ncoe;                   % loop over the coefficient
  for ilay  = 1:nlay;
    clear coefX cXvq vq;
    coefX = coef(:,ilay,ic);      % figure(1);clf;plot(fchan,coefX,'.-');grid on;

    % regularize the grid of the coefficient 
    vq = interp1(fchan,coefX,fiasi);
    vq(isnan(vq)) = 0.0;
       % whos vq fchan coefX fiasi
    
    % convert the coefficient from the regularized IASI grid to CrIS
    [cXvq, cfreq] = iasi2cris(vq, fiasi, opt1);
       % whos cXvq cfreq fcris

    % record common channels
    [si sj] = seq_match(cfreq, fchan);
  
    ccoef(:,ilay, ic) = real(cXvq(si));
  end
  fprintf(1,'.');
end
fprintf(1,'\n');
   whos si sj cfreq ccoef;

%{
ic = 1; ilay = 90;
figure(1);clf;plot(fchan,coef(:,ilay,ic),'.-',cfreq(si),ccoef(:,ilay,ic),'.-');grid on;

figure(4);clf;
  h1=subplot(2,1,1);plot(fchan,coefX,'.-',cfreq,real(cXvq),'.-');grid on;xlim([640 800]);
  ylabel('coef arb.');title('iasi2cris set1 lay20 c8');
  h2=subplot(2,1,2);plot(fchan(sj), (coefX(sj) - real(cXvq(si)))./coefX(sj),'.-'); 
  grid on; xlim([640 800]);ylabel('fract diff');xlabel('wn cm-1');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
  set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
  pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
%}

% ---------------------------------------------
% Preparation to write binary coefficient file
% ---------------------------------------------
o_info = struct;
o_info.set       = cset;
o_info.ncoef     = ncoe;
o_info.nlay      = info.nlay;
o_info.gasid     = info.gasid;
o_info.nbreakout = info.nbreakout;
o_info.ncoefeach = info.ncoefeach;
o_info.startin   = info.startind;
%%o_info.lcon      = info.lcon;
o_info.nchan     = numel(si);
if(mset == 9)                   % extra stuff for optran
o_info.navgpred  = info.navgpred;
o_info.avgpred   = info.avgpred;
o_info.ogrid     = info.ogrid;

% Open output file
opath   = '/asl/s1/chepplew/data/sarta_database/Data_i2c_Mar16/Coef/';
clear ofiln;
if(mset >= 1 && mset <= 7 || mset == 9) 
  ofiln   = ['i2c_set' sprintf('%d',mset) '.dat'];
elseif (mset == 8)
  ofiln   = ['i2c_set' sprintf('%d',mset) '_gas' sprintf('%d',info.gasid) '.dat'];
elseif (mset == 10)
  ofiln   = ['i2c_set' sprintf('%d',mset) '_gas' sprintf('%d',info.gasid) '.dat'];  
end
disp([num2str(mset) '  ' ofiln]);

fid = fopen(strcat(opath, ofiln),'w','ieee-be');

% Write merged OPTRAN header
if (mset == 9 & isfield(o_info,'navgpred'))
   ifm = 4*o_info.nlay;
   icount=fwrite(fid,ifm,'integer*4');
   icount=fwrite(fid,o_info.ogrid,'real*4');
   icount=fwrite(fid,ifm,'integer*4');
   for ip=1:o_info.navgpred
      icount=fwrite(fid,ifm,'integer*4');
      icount=fwrite(fid,o_info.avgpred(:,ip),'real*4');
      icount=fwrite(fid,ifm,'integer*4');
   end
end


% Value of FORTRAN record marker for each channel
ifm = round( 4*(1 + 1 + o_info.nlay*o_info.ncoef) ); % exact integer

% Loop over the channels
for k = 1:o_info.nchan

   % FORTRAN start-of-record marker
   icount = fwrite(fid,ifm,'integer*4');

   % data for this channel
   icount = fwrite(fid,si(k),'integer*4');
   icount = fwrite(fid,cfreq(si(k)),'real*4');
   for i = 1:o_info.nlay
      icount = fwrite(fid,ccoef(k,i,:),'real*4');
   end

   % FORTRAN end-of-record marker
   icount=fwrite(fid,ifm,'integer*4');

end % end of loop over bands

fclose(fid);

iok = info.nchan;

%{
%zconts = 1.e-7*[-3:0.2:0.4];
%figure(5);clf;contourf([fchan,coef(:,:,1)],zconts);
%for sets 1:8, nlay = 100. set 9: nlay = 300.[XY YY] = meshgrid(fchan,1:300);
[XY YY] = meshgrid(fchan, 1:300);
[XZ YZ] = meshgrid(cfreq(si), 1:300);
for ic = 1:ncoe;
figure(5); clf; subplot(2,1,1);H1 = surf(XY, YY, coef(:,:,ic)');colorbar; set(H1,'edgecolor','none');
title(['IASI coef gas ' num2str(info.gasid) ' set ' num2str(cset) ' num ' num2str(ic)]);
subplot(2,1,2);H2 = surf(XZ, YZ, real(ccoef(:,:,ic)'));colorbar; set(H2,'edgecolor','none');
  title('iasi 2 cris');xlabel('freq');ylabel('layer');
    %saveas(gcf,['./figs/i2c_set' num2str(cset) '_gas' num2str(info.gasid) '_coef' num2str(ic) '_surf.png'],'png');
pause(5);
end

% for set 11, ntc_7term.dat, there is only 1 layer, 7 coeffs. therm.dat has 6 coeffs
[XY YY] = meshgrid(fchan,1:6);
[XZ YZ] = meshgrid(cfreq(si), 1:6);
clear junk; junk = squeeze(coef(:,1,:));
figure(5);clf;subplot(2,1,1);H1=surf(XY, YY, junk');colorbar; set(H1,'edgecolor','none');
  title('IASI coef therm wn vs coef');
clear junk; junk = squeeze(ccoef(:,1,:));
subplot(2,1,2);H2=surf(XZ, YZ, junk');colorbar; set(H2,'edgecolor','none');
  title('iasi2cris');xlabel('freq');ylabel('num coef');
  %saveas(gcf,'./figs/i2c_set11_therm_coef_surf.png','png');
%}
