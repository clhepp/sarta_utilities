function [ichan, fchan, coef]=rd_concoef( fname );

% function [ichan, fchan, coef]=rd_concoef( fname );
%
% Read in fast model water continuum coefficients
%
% Input:
%    fname : {string} name of water continuum coefficient file
%
% Output:
%    ichan : [1 x nchan] channel ID number
%    fchan : [1 x nchan] channel frequency
%    coef : [ncoef x nlay x nchan] coefficients
%

% Created: Scott Hannon, 29 May 2002 - based on interp_watercon.m
% Update: 05 Jan 2005, S.Hannon - pre-declare arrays; more error checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncoef = 7;
nlay = 100;


% Check file exists
d = dir(fname);
if (length(d) ~= 1)
   error(['Problem opening file: ' fname]);
end
nbytes=d.bytes;

% Open binary file
% Note: FORTRAN format big-endian IEEE 4-byte reals
fid=fopen(fname,'r','ieee-be');


% Determine number of channels from expected record size and
% actual file size
% Note: FORTRAN format big-endian IEEE 4-byte reals
xirecsize=4*(2 + nlay*ncoef); % Expected record size
nchan=nbytes/(xirecsize + 8);  % plus 8=4*2 for the record markers
if (nchan ~= round(nchan))
   nchan
   error('Unexpected total file size based on expected record size')
end


% Declare output arrays
ichan = zeros(1,nchan);
fchan = zeros(1,nchan);
coef = zeros(ncoef,nlay,nchan);


% Loop over channels
for ii=1:nchan

   irecsize=fread(fid,1,'integer*4');
   if (irecsize ~= xirecsize)
     irecsize
     xirecsize
     error('Unexpected record size')
   end

   ichan(ii)=fread(fid,1,'integer*4');
   fchan(ii)=fread(fid,1,'real*4');

   for il=1:nlay
      coef(:,il,ii)=fread(fid,ncoef,'real*4');
   end
   irecsize=fread(fid,1,'integer*4');

   if (mod(ii,100)== 0)
     disp(['doing channel ' int2str(ii)])
   end

end

ii=fclose(fid);

%%% end of program %%%
