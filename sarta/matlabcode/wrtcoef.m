function [iok] = wrtcoef(ichan, fchan, coef, info, fname);

% function [iok] = wrtcoef(ichan, fchan, coef, info, fname);
%
% Write a binary FORTRAN coefficient data file (set1-7). Also
% see "rdcoef.m".
%
% Input:
%    ichan   = [nchan x 1] channel ID (JPL numbering)
%    fchan   = [nchan x 1] center frequency (wavenumber)
%    coef    = [nchan x nlay x ncoef] coeffcients
%    info    = structure of supplemental info
%    fname   = {string} name of binary FORTRAN data file to write
%
% Output:
%    iok     = [1 x 1] 0 if error, otherwise nchan if OK 
%

% Created: 02 January 2008, Scott Hannon - based on rdcoef.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iok=0;

if (isfield(info,'set'))
   set=info.set;
else
   error('input structure "info" lacks field "set"')
end

% Open output file
fid=fopen(fname,'w','ieee-be');


%%%
% Write merged OPTRAN header
if (set == 9 & isfield(info,'navgpred'))
   ifm = 4*info.nlay;
   icount=fwrite(fid,ifm,'integer*4');
   icount=fwrite(fid,info.ogrid,'real*4');
   icount=fwrite(fid,ifm,'integer*4');
   for ip=1:info.navgpred
      icount=fwrite(fid,ifm,'integer*4');
      icount=fwrite(fid,info.avgpred(:,ip),'real*4');
      icount=fwrite(fid,ifm,'integer*4');
   end
end
%%%

% Value of FORTRAN record marker for each channel
ifm=round( 4*(1 + 1 + info.nlay*info.ncoef) ); % exact integer


% Loop over the channels
for ic=1:info.nchan

   % FORTRAN start-of-record marker
   icount=fwrite(fid,ifm,'integer*4');

   % data for this channel
   icount=fwrite(fid,ichan(ic),'integer*4');
   icount=fwrite(fid,fchan(ic),'real*4');
   for il=1:info.nlay
      icount=fwrite(fid,coef(ic,il,:),'real*4');
   end

   % FORTRAN end-of-record marker
   icount=fwrite(fid,ifm,'integer*4');

end % end of loop over bands

fclose(fid);

iok = info.nchan;

%%% end of function %%%
